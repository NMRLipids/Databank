#!/usr/bin/env python
# coding: utf-8

import os
import sys
import urllib.request

import urllib.request
from urllib.error import URLError,HTTPError,ContentTooShortError

import socket

sys.path.append('..')
from DatabankLib.databankLibrary import lipids_dict, databank, loadMappingFile
from DatabankLib.databankio import resolve_download_file_url
import DatabankLib.form_factor as form_factor

path = '../../Data/Simulations/'

db_data = databank(path)
systems = db_data.get_systems()
 
for system in systems:
    software = system['SOFTWARE']
    # download trajectory and gro files
    system_path = '../../Data/Simulations/' +  system['path']
    doi = system.get('DOI')
    skipDownloading: bool = (doi == 'localhost')
    if skipDownloading:
        print("NOTE: The system with 'localhost' DOI should be downloaded by the user.")

    if (os.path.isfile(system_path + "/FormFactor.json")):
        continue

    if system['TYPEOFSYSTEM'] == 'miscellaneous':
        continue
    
    print('Analyzing', system_path)

    try:
        if system['UNITEDATOM_DICT']:
            #continue
            print('United atom simulation')
    except:
        pass

    try:
        if system['WARNINGS']['ORIENTATION']:
            print('Skipping due to ORIENTATION warning:', system['WARNINGS']['ORIENTATION'])#, system['WARNINGS']['PBC'])
            continue
    except:
        pass

    try:
        if system['WARNINGS']['PBC'] == 'hexagonal-box':
            print('Skipping due to PBC warning:', system['WARNINGS']['PBC'])#, system['WARNINGS']['PBC'])
            continue
    except:
        pass

    try:
        if system['WARNINGS']['NOWATER']:
            print('Skipping because there is not water in the trajectory.')
            continue
    except:
        pass
    
    output_name = ""
    
    trj_name = system_path + system['TRJ'][0][0]

    socket.setdefaulttimeout(15)

    if (skipDownloading):
        if (not os.path.isfile(trj_name)):
            raise FileNotFoundError(f"Trajectory should be downloaded [{trj_name}] by user")
    else:
        trj_url = resolve_download_file_url(system['DOI'], system['TRJ'][0][0])
        if (not os.path.isfile(trj_name)):
            print('Downloading trajectory with the size of ', system['TRAJECTORY_SIZE'], ' to ', system['path'])
            response = urllib.request.urlretrieve(trj_url, trj_name)

    # make a function like this
    if 'gromacs' in software:
        tpr_name = system_path + system['TPR'][0][0]
        
        if (skipDownloading):
            if (not os.path.isfile(tpr_name)):
                raise FileNotFoundError(f"TPR should be downloaded [{tpr_name}] by user")
        else:
            tpr_url = resolve_download_file_url(doi, system['TPR'][0][0])
            if (not os.path.isfile(tpr_name)):
                response = urllib.request.urlretrieve(tpr_url, tpr_name)
               
    if 'openMM' in software or 'NAMD' in software:
        pdb = system.get('PDB')
        pdb_name = system_path + pdb[0][0]
        if (skipDownloading):
            if (not os.path.isfile(pdb_name)):
                raise FileNotFoundError(f"PDB should be downloaded [{pdb_name}] by user")
        else:
            pdb_url = resolve_download_file_url(doi, pdb[0][0])
            if (not os.path.isfile(pdb_name)):
                response = urllib.request.urlretrieve(pdb_url, pdb_name)
                        
    EQtime = float(system['TIMELEFTOUT'])*1000

    # FIND LAST CARBON OF SN-1 TAIL AND G3 CARBON
    for molecule in system['COMPOSITION']:
        if molecule in lipids_dict:
            mapping_file = system['COMPOSITION'][molecule]['MAPPING']
            mapping = loadMappingFile(mapping_file)

            #TODO: rewrite via lipid dictionary!
            for nm in ["M_G3_M", "M_G13_M", "M_C32_M"]:
                try:
                    G3atom = mapping[nm]['ATOMNAME']
                    continue
                except:
                    pass
            
            #TODO: rewrite via lipid dictionary
            if "M_G1C4_M" in mapping.keys():
                for Cindex in range(4,30):
                    atom = 'M_G1C' + str(Cindex) + '_M'
                    try:
                        lastAtom = mapping[atom]['ATOMNAME']
                    except:
                        continue
            elif "M_G11C4_M" in mapping.keys():
                for Cindex in range(4,30):
                    atom = 'M_G11C' + str(Cindex) + '_M'
                    try:
                        lastAtom = mapping[atom]['ATOMNAME']
                    except:
                        continue
            elif "M_CA4_M" in mapping.keys():
                for Cindex in range(4,30):
                    atom = 'M_CA' + str(Cindex) + '_M'
                    try:
                        lastAtom = mapping[atom]['ATOMNAME']
                    except:
                        continue
                
    print(lastAtom, G3atom)

    # Center around one lipid tail CH3 to guarantee all lipids in the same box
    if 'gromacs' in system['SOFTWARE']:

        if ('WARNINGS' in system and
            system['WARNINGS'] is not None and
            'GROMACS_VERSION' in system['WARNINGS'] and
            system['WARNINGS']['GROMACS_VERSION'] == 'gromacs3'):
            trjconvCOMMAND = 'trjconv'
            makendxCOMMAND = 'make_ndx'
        else:
            trjconvCOMMAND = 'gmx trjconv'
            makendxCOMMAND = 'gmx make_ndx'
        
        os.system('rm foo.ndx')
        os.system(f'echo "a {lastAtom}\nq" | {makendxCOMMAND} -f {tpr_name} -o foo.ndx')
        os.system("tail -n1 foo.ndx | awk '{print $NF}'")
        os.system('echo "[ centralAtom ]" >> foo.ndx')
        os.system("tail -n2 foo.ndx | head -n1 |  awk '{print $NF}' >> foo.ndx")

        xtcwhole= system_path + '/whole.xtc'
        xtcfoo = system_path + '/foo2.xtc'
        xtccentered= system_path + '/centered.xtc'
        if (not os.path.isfile(xtccentered)):
            print("Make molecules whole in the trajectory")
            #if unitedAtom and system['TRAJECTORY_SIZE'] > 15000000000:
            #    print("United atom trajectry larger than 15 Gb. Using only every third frame to reduce memory usage.")
            #    os.system('echo System | gmx trjconv -f ' + trj_name + ' -s ' + tpr_name + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime) + ' -skip 3')
            #else:
            if (not os.path.isfile(xtcwhole)):
                os.system(f'echo System |  {trjconvCOMMAND} -f {trj_name} -s {tpr_name} -o {xtcwhole} -pbc mol -b {str(EQtime)}')

            if (not os.path.isfile(xtcfoo)):
                os.system(f'echo "centralAtom\nSystem" |  {trjconvCOMMAND} -center -pbc mol -n foo.ndx -f {xtcwhole} -s {tpr_name} -o {xtcfoo}') 

            os.system('rm foo.ndx')
            os.system('rm ' + xtcwhole)
                
        #Center around the center of mass of all the g_3 carbons
        #if (not os.path.isfile(xtccentered)):        
            os.system(f'echo "a {G3atom}\nq" | {makendxCOMMAND} -f {tpr_name} -o foo.ndx')
            os.system(f'echo "{G3atom}\nSystem" |  {trjconvCOMMAND} -center -pbc mol -n foo.ndx -f {xtcfoo} -s {tpr_name} -o {xtccentered}')
            os.system('rm ' + xtcfoo)            
    else:
        print('Centering for other than Gromacs may not work if there are jumps over periodic boundary conditions in z-direction.')

    if (not os.path.isfile(system_path + "/FormFactor.json")):
        try:
            if 'gromacs' in system['SOFTWARE']:
                form_factor.FormFactor(system_path, tpr_name, xtccentered, 200, output_name,  system)
            if 'openMM' in system['SOFTWARE'] or 'NAMD' in system['SOFTWARE']:
                form_factor.FormFactor(system_path, pdb_name, trj_name, 200, output_name,  system)
        except ValueError as e:
            # Here it was expected to have allow_pickle-type errors.
            # but I suppose, we cannot simply ignore them because it means that the code breaks at this place,
            # and the running user should be informed about it!
            raise e