import os, sys
import json, yaml
import urllib.request
import re
import buildh

sys.path.append('..')
from DatabankLib.databankLibrary import *
from DatabankLib.databankio import resolve_download_file_url

### Initializing the databank
databankPath = '../../'
systems = initialize_databank(databankPath)

## intialize counters for analyzed systems
ready = 0
skipped = 0

## loop over simulations
for system in systems:
    path = system['path']

    ## Check if order parameters are calculated or something in the system prevents order parameter calculation
    for key in system['COMPOSITION']:
        outfilename = databankPath + '/Data/Simulations/' + path + key + 'OrderParameters.json'
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
            continue

    if FileFound:
        skipped += 1
        continue

    ## If order parameters are not calculated and system ok, then continue to calculating order parameters
    print('Analyzing: ', path)


    ## Download the simulations and create MDAnalysis universe (the universe is not used here but this script downloads the data)
    u = system2MDanalysisUniverse(system)

    ## Store the local path of the trajectory
    trj_name = databankPath + '/Data/Simulations/' + path + system.get('TRJ')[0][0]

    ## Software and time for equilibration period
    software=system['SOFTWARE']
    EQtime=float(system['TIMELEFTOUT'])*1000

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

    ## Set topology file names and make Gromacs trajectories whole
    if 'gromacs' in software:
        tpr_name = databankPath + '/Data/Simulations/' + path + system.get('TPR')[0][0]

        xtcwhole = databankPath + '/Data/Simulations/' + path + '/whole.xtc'
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
        pdb_name = databankPath + '/Data/Simulations/' + path + system.get('PDB')[0][0]
        if (not os.path.isfile(pdb_name)):
            pdb_url = resolve_download_file_url(system.get('DOI'), pdb_name)
            response = urllib.request.urlretrieve(pdb_url, pdb_name)
    else:
        print('Order parameter calculation for other than gromacs, openMM and NAMD are yet to be implemented.')
        continue

    ## Calculate order parameters
    if unitedAtom and 'gromacs' in software:
        topfile = databankPath + '/Data/Simulations/' + path + '/frame0.gro'
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
            
            def_fileNAME = databankPath + '/Data/Simulations/' + path + key + '.def' 
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
            ordPfile = databankPath + '/Data/Simulations/' + path + key + 'OrderParameters.dat' 

            lipid_json_file = ['./lipid_json_buildH/' + system['UNITEDATOM_DICT'][key] + '.json']

            if (not os.path.isfile(lipid_json_file[0])):
                lipid_json_file = None
            
            print(system['UNITEDATOM_DICT'][key])
            buildh.launch(coord_file=topfile, def_file=def_fileNAME, lipid_type=system['UNITEDATOM_DICT'][key], lipid_jsons=lipid_json_file, traj_file=xtcwhole , out_file=f"{ordPfile}.buildH", ignore_CH3s=True)

            outfile = open(ordPfile,'w')
            outfile.write("Atom     Average OP     OP stem\n")
        
            data = {}
            outfile2 = databankPath + '/Data/Simulations/' + path + key + 'OrderParameters.json'
        
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
            gro = databankPath + '/Data/Simulations/' + path + '/conf.gro'
            
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
                outfilename = databankPath + '/Data/Simulations/' + path + key + 'OrderParameters.dat'
                outfilename2 = databankPath + '/Data/Simulations/' + path + key + 'OrderParameters.json'
                if (os.path.isfile(outfilename2)):
                    print('Order parameter file already found')
                    continue
                outfile = open(outfilename,'w')

                if 'gromacs' in software:
                    try:
                        OrdParam = find_OP(mapping_file,tpr_name,xtcwhole,resname)
                    except:
                        print('Using tpr did not work, trying with gro')
                        OrdParam = find_OP(mapping_file,gro,xtcwhole,resname)

                if 'openMM' in software or 'NAMD' in software:
                    OrdParam = find_OP(mapping_file,pdb_name,trj_name,resname)
                        
                outfile.write("Atom     Average OP     OP stem\n")
    
                data = {}
                outfile2 = outfilename2 

                for i,op in enumerate(OrdParam):
                    resops = op.get_op_res
                    (op.avg, op.std, op.stem) = op.get_avg_std_stem_OP
                    outfile.write(f'{op.name} {str(op.avg)} {str(op.stem)}\n')
    
                    data[str(op.name)]=[]
                    data[str(op.name)].append(op.get_avg_std_stem_OP)
        
                with open(outfile2, 'w') as f:
                    json.dump(data,f)
                outfile.close()
                f.close()
    
        print("Order parameters calculated and saved to ",path)

    ready = ready + 1
        
print('Order parameters calculated for ', ready , 'systems.')
print('Already calculated order parameters found for', skipped , 'systems.')
