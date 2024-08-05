
import json
import os
import logging
import re
import traceback
import buildh
import urllib.request
from tqdm import tqdm

logger = logging.getLogger(__name__)

from . import *
from .databank_defs import *
from .core import *
from .databankLibrary import GetNlipids, loadMappingFile, system2MDanalysisUniverse
from .jsonEncoders import CompactJSONEncoder
from .databankio import resolve_download_file_url
from .databankop import find_OP

## Universal success codes

RCODE_SKIPPED: int = 0
""" Success code 0: calculation skipped """

RCODE_COMPUTED: int = 1
""" Success code 1: calculation sucessful """

RCODE_ERROR: int = 2
""" Success code 2: calculation failed """


#TODO: use in calcAPL.py
def computeAPL(system: dict, recompute: bool = False) -> int:
    """Generate apl.json analysis file for a system.

    Args:
        system (dict): one of systems of the Databank
        recompute (bool, optional): Delete previous apl.json and recompute it if True. Defaults to False.
    Returns:
        int success code (RCODE_...)
    """
    ## reading software and file path for simulations
    software = system['SOFTWARE']
    path = system['path']

    ## this is checking if area per lipid is already calculated for the systems
    outfilename = os.path.join(NMLDB_SIMU_PATH,  path, 'apl.json')
    if os.path.isfile(outfilename):
        if recompute:
            os.unlink(outfilename)
        else:
            return RCODE_SKIPPED
    
    print('Analyzing: ', path)
    print('Will write into: ', outfilename)

    ## calculates the total number of lipids
    Nlipid = GetNlipids(system)

    ## makes MDAnalysis universe from the system. This also downloads the data if not yet locally available
    u = system2MDanalysisUniverse(system)

    if u is None:
        print('Generation of MDAnalysis universe failed in folder', path)
        return RCODE_ERROR
    
    ## this calculates the area per lipid as a function of time and stores it in the databank
    apl = {}
    for ts in tqdm(u.trajectory, desc='Scanning the trajectory'):
        if u.trajectory.time >= system['TIMELEFTOUT']*1000:
            dimensions = u.dimensions
            aplFrame = u.dimensions[0]*u.dimensions[1]*2/Nlipid
            apl[u.trajectory.time] = aplFrame

    with open(outfilename, 'w') as f:
        json.dump(apl, f, cls=CompactJSONEncoder)
    
    return RCODE_COMPUTED

#TODO: implement onlyLipid
def computeOP(system: dict, recompute: bool = False) -> int:
    """_summary_

    Args:
        system (dict): _description_
        recompute (bool, optional): _description_. Defaults to False.

    Returns:
        int: _description_
    """
    path = system['path']

    ## Check if order parameters are calculated or something in the system prevents order parameter calculation
    for key in system['COMPOSITION']:
        outfilename = os.path.join(NMLDB_SIMU_PATH, path, key + 'OrderParameters.json')
        #print(outfilename)
        if ( os.path.isfile(outfilename)  or 
             ('WARNINGS' in system.keys() and 
              type(system['WARNINGS']) is dict and
              'AMBIGUOUS_ATOMNAMES' in system['WARNINGS'].keys() and 
              key in system['WARNINGS']['AMBIGUOUS_ATOMNAMES'])
           ):
            FileFound = True
        elif key in lipids_dict:
            FileFound = False
            break

    if FileFound and not recompute:
        return RCODE_SKIPPED

    ## If order parameters are not calculated and system ok, then continue to calculating order parameters
    print('Analyzing: ', path)

    ## Download the simulations and create MDAnalysis universe (the universe is not used here but this script downloads the data)
    u = system2MDanalysisUniverse(system)

    ## Store the local path of the trajectory
    trj_name = os.path.join(NMLDB_SIMU_PATH, path, system.get('TRJ')[0][0])

    ## Software and time for equilibration period
    software = system['SOFTWARE']
    EQtime = float(system['TIMELEFTOUT'])*1000

    ## Check if all or united atom simulation
    try:
        unitedAtom = system['UNITEDATOM_DICT']
    except:
        unitedAtom = False

    ## Check relevant warnings
    g3switch = ( 'WARNINGS' in system and 
                 type(system['WARNINGS']) is dict and
                 'GROMACS_VERSION' in system['WARNINGS'] and 
                 system['WARNINGS']['GROMACS_VERSION'] == 'gromacs3' )
    trjconvCOMMAND = 'trjconv' if g3switch else 'gmx trjconv'

    try:
        ## Set topology file names and make Gromacs trajectories whole
        if 'gromacs' in software:
            tpr_name = os.path.join(NMLDB_SIMU_PATH, path, system.get('TPR')[0][0])

            xtcwhole = os.path.join(NMLDB_SIMU_PATH, path, 'whole.xtc')
            if (not os.path.isfile(xtcwhole)):
                execStr = f'echo System | {trjconvCOMMAND} -f {trj_name} -s {tpr_name} -o {xtcwhole} -pbc mol -b {str(EQtime)}'
                print("Make molecules whole in the trajectory")
                if unitedAtom and system['TRAJECTORY_SIZE'] > 15e9:
                    print("United atom trajectry larger than 15 Gb. Using only every third frame to reduce memory usage.")
                    execStr += " -skip 3"
                rCode = os.system(execStr)
                if (rCode != 0):
                    raise RuntimeError("trjconv exited with error (see above)")
        elif 'openMM' in software or 'NAMD' in software:
            pdb_name = os.path.join(NMLDB_SIMU_PATH, path, system.get('PDB')[0][0])
            if (not os.path.isfile(pdb_name)):
                pdb_url = resolve_download_file_url(system.get('DOI'), pdb_name)
                response = urllib.request.urlretrieve(pdb_url, pdb_name)
        else:
            print('Order parameter calculation for other than gromacs, openMM and NAMD are yet to be implemented.')
            return RCODE_ERROR

        ## Calculate order parameters
        if unitedAtom and 'gromacs' in software:
            topfile = os.path.join(NMLDB_SIMU_PATH, path, 'frame0.gro')
            if g3switch:
                rCode = os.system(f'echo System | editconf -f {tpr_name} -o {topfile}')
                if (rCode != 0):
                    raise RuntimeError("editconf exited with error (see above)")
            else:
                rCode = os.system('echo System | {trjconvCOMMAND} -f {xtcwhole} -s {tpr_name} -dump 0 -o {topfile}' )
                if (rCode != 0):
                    raise RuntimeError(f"trjconv ({trjconvCOMMAND}) exited with error (see above)") 

            for key in system['UNITEDATOM_DICT']:
            #construct order parameter definition file for CH bonds from mapping file
                mapping_file = system['COMPOSITION'][key]['MAPPING']
                mapping_dict = loadMappingFile(mapping_file)
                
                def_fileNAME = os.path.join(NMLDB_SIMU_PATH, path, key + '.def')
                def_file = open(def_fileNAME, 'w')
                
                previous_line = ""            

                regexp1_H = re.compile(r'M_[A-Z0-9]*C[0-9]*H[0-9]*_M')
                regexp2_H = re.compile(r'M_G[0-9]*H[0-9]*_M')
                regexp1_C = re.compile(r'M_[A-Z0-9]*C[0-9]*_M')
                regexp2_C = re.compile(r'M_G[0-9]_M')
                
                for mapping_key in mapping_dict.keys():
                    if regexp1_C.search(mapping_key) or regexp2_C.search(mapping_key):
                        atomC = [mapping_key, mapping_dict[mapping_key]['ATOMNAME']]
                        atomH = []
                    elif regexp1_H.search(mapping_key) or regexp2_H.search(mapping_key):
                        atomH = [mapping_key, mapping_dict[mapping_key]['ATOMNAME']]
                    else:
                        atomC = []
                        atomH = []

                    if atomH:
                        items = [atomC[1], atomH[1], atomC[0], atomH[0]]
                        def_line = items[2] + "&" + items[3] + " " + key + " " + items[0] + " " + items[1] + "\n"
                        if def_line != previous_line:
                            def_file.write(def_line)
                            previous_line = def_line
                def_file.close()            
                
                #Add hydrogens to trajectory and calculate order parameters with buildH
                ordPfile = os.path.join(NMLDB_SIMU_PATH, path, key + 'OrderParameters.dat')

                lipid_json_file = ['./lipid_json_buildH/' + system['UNITEDATOM_DICT'][key] + '.json']

                if (not os.path.isfile(lipid_json_file[0])):
                    lipid_json_file = None
                
                print(system['UNITEDATOM_DICT'][key])
                buildh.launch(coord_file=topfile, def_file=def_fileNAME, lipid_type=system['UNITEDATOM_DICT'][key], lipid_jsons=lipid_json_file, traj_file=xtcwhole , out_file=f"{ordPfile}.buildH", ignore_CH3s=True)

                outfile = open(ordPfile,'w')
                outfile.write("Atom     Average OP     OP stem\n")
            
                data = {}
                outfile2 = os.path.join(NMLDB_SIMU_PATH, path, key + 'OrderParameters.json')
            
                with open(ordPfile + '.buildH') as OPfile:
                    lines = OPfile.readlines()
                    for line in lines:
                        if "#" in line:
                            continue
                        line2 = line.split()[0].replace('&',' ') + "  " + line.split()[4] + "  " + line.split()[5] + " " + line.split()[6] + "\n"
                        outfile.write(line2)

                        OPname = line.split()[0].replace('&',' ') #line.split()[0] + " " + line.split()[1]
                        OPvalues = [line.split()[4], line.split()[5] ,line.split()[6]]
                        data[str(OPname)]=[]
                        data[str(OPname)].append(OPvalues)
            
                with open(outfile2, 'w') as f:
                    json.dump(data,f)

                outfile.close()

        else:
            if 'gromacs' in software:
                gro = os.path.join(NMLDB_SIMU_PATH, path, 'conf.gro')
                
                #make gro file
                print("\n Makin gro file")
                if g3switch:
                    rCode = os.system(f'echo System | editconf -f {tpr_name} -o {gro}') 
                    if (rCode != 0):
                        raise RuntimeError("editconf exited with error (see above)")
                else:
                    rCode = os.system(f'echo System | gmx trjconv -f {trj_name} -s {tpr_name} -dump 0 -o {gro}')
                    if (rCode != 0):
                        raise RuntimeError("trjconv exited with error (see above)")
                        
            for key in system['COMPOSITION']:

                if ( 'WARNINGS' in system.keys() and 
                    system['WARNINGS'] is not None and
                    'AMBIGUOUS_ATOMNAMES' in system['WARNINGS'].keys() and 
                    key in system['WARNINGS']['AMBIGUOUS_ATOMNAMES'] ):
                    print(path, key)
                    print('Order parameters cannot be calculated if atom names are ambiguous.')
                    continue
                
                if key in lipids_dict.keys():
                    print('Calculating ', key,' order parameters')
                    mapping_file = system['COMPOSITION'][key]['MAPPING']
                    resname = system['COMPOSITION'][key]['NAME']
                    outfilename = os.path.join(NMLDB_SIMU_PATH, path, key + 'OrderParameters.dat')
                    outfilename2 = os.path.join(NMLDB_SIMU_PATH, path, key + 'OrderParameters.json')
                    if (os.path.isfile(outfilename2)):
                        print('Order parameter file already found')
                        continue

                    if 'gromacs' in software:
                        try:
                            OrdParam = find_OP(mapping_file, tpr_name, xtcwhole, resname)
                        except Exception as e:
                            logger.warning("We got this exception: \n    " + str(e) )
                            logger.warning('But we will try rebuild the Universe from GROM if using tpr did not work!')
                            OrdParam = find_OP(mapping_file,gro,xtcwhole,resname)

                    if 'openMM' in software or 'NAMD' in software:
                        OrdParam = find_OP(mapping_file,pdb_name,trj_name,resname)
                    
                    data = {}

                    with open(outfilename,'w') as outfile:
                        outfile.write("Atom     Average OP     OP stem\n")

                        for i,op in enumerate(OrdParam):
                            resops = op.get_op_res
                            (op.avg, op.std, op.stem) = op.get_avg_std_stem_OP
                            outfile.write(f'{op.name} {str(op.avg)} {str(op.stem)}\n')
            
                            data[str(op.name)]=[]
                            data[str(op.name)].append(op.get_avg_std_stem_OP)
            
                    with open(outfilename2, 'w') as f:
                        json.dump(data, f, cls=CompactJSONEncoder)
        
        print("Order parameters calculated and saved to ", path)

    except Exception as e:
        print('Calculation failed for ' +  system['path'])
        print(str(e))
        print(traceback.format_exc())
        return RCODE_ERROR

    return RCODE_COMPUTED