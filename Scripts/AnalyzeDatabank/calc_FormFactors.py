import os
import sys
import numpy as np
import json
import matplotlib.pyplot as plt
import MDAnalysis
import urllib.request
import yaml

import urllib.request
from urllib.error import URLError,HTTPError,ContentTooShortError

# From time monitoring
from tqdm import tqdm

import socket

sys.path.insert(1, '../BuildDatabank/')
from databankLibrary import download_link, lipids_dict, databank, read_mapping_file
import form_factor

path = '../../Data/Simulations/'

db_data = databank(path)
systems = db_data.get_systems()

with open("system_already_run_17_11_21.out","w") as f:
    f.write("Listing analysed systems \n \n ")   
    
     
for system in systems:
    software=system['SOFTWARE']
    #download trajectory and gro files
    system_path = '../../Data/Simulations/' +  system['path']
    doi = system.get('DOI')

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

    
    #try:
    #    if system['WARNINGS']['PBC']:
    #        print('Skipping due to PBC warning')#, system['WARNINGS']['PBC'])
    #        continue
    #except:
    #    pass

    
    output_name = ""
    
    trj_name = system_path + system['TRJ'][0][0]
    trj_url = download_link(system['DOI'], system['TRJ'][0][0])

    socket.setdefaulttimeout(15)

    download_failed = False

    
    # make a function like this
    if 'gromacs' in software:
        tpr_name = system_path + system['TPR'][0][0]
        tpr_url = download_link(system['DOI'], system['TPR'][0][0])

        if (not os.path.isfile(tpr_name)):
            try:
                url_size = urllib.request.urlopen(tpr_url).length
                response = urllib.request.urlretrieve(tpr_url, tpr_name)
            
            except ContentTooShortError:
                download_failed = True
                print("Content too short error.")
            except HTTPError as e:
                download_failed = True
                print(e)
            except URLError as ue:
                download_failed = True
                print("failed to download")
            except socket.timeout as se:
                download_failed = True
                print("socket time out")
            except Exception as ee:
                download_failed = True
                print(ee)
                #check if the file is fully downloaded
            else:
                size = os.path.getsize(tpr_name)
                print("size of the file "  + " to be downloaded: " + str(url_size))
                print("size of the file " + " after download: " + str(size) )
                if url_size != size:
                    print("Download of the file "  + " was interrupted.")
               
    if 'openMM' in software or 'NAMD' in software:
        #print(system)
        #print(software)
        pdb = system.get('PDB')
        pdb_name = system_path + system.get('PDB')[0][0]
        pdb_url = download_link(doi, pdb[0][0])
        if (not os.path.isfile(pdb_name)):
            response = urllib.request.urlretrieve(pdb_url, pdb_name)


                        
    if (not os.path.isfile(trj_name)) and download_failed == False:
        
        try:
            url_size = urllib.request.urlopen(trj_url).length
            response = urllib.request.urlretrieve(trj_url, trj_name)
       
        except ContentTooShortError:
            download_failed = True
            print("Content too short error.")
        except HTTPError as e:
            download_failed = True
            print(e)
        except URLError as ue:
            download_failed = True
            print("failed to download")
        except socket.timeout as se:
            download_failed = True
            print("socket time out")
        except Exception as ee:
            download_failed = True
            print(ee)
        else:
            #check if the file is fully downloaded
            size = os.path.getsize(trj_name)
            print("size of the file " +  " to be downloaded: " + str(url_size))
            print("size of the file " +  " after download: " + str(size) )
            if url_size != size:
                print("Download of the file "  + " was interrupted.")
            print("Download was a success!!!!!!")
     
    if download_failed == True:
        with open("system_already_run_17_11_21.out","a") as f:
            f.write(system_path+" Download failed !!!\n")    


    EQtime=float(system['TIMELEFTOUT'])*1000

    # FIND LAST CARBON OF SN-1 TAIL AND G3 CARBON
    for molecule in system['COMPOSITION']:
        if molecule in lipids_dict:
            mapping_file = '../BuildDatabank/mapping_files/' + system['COMPOSITION'][molecule]['MAPPING']
            with open(mapping_file, "r") as yaml_file:
                mapping = yaml.load(yaml_file,  Loader=yaml.FullLoader)
            try:
                G3atom = mapping['M_G3_M']['ATOMNAME']
            except:
                try:
                    G3atom = mapping['M_C32_M']['ATOMNAME']
                except:
                    pass
        
            for Cindex in range(1,30):
                atom = 'M_G1C' + str(Cindex) + '_M'
                try:
                    lastAtom = mapping[atom]['ATOMNAME']
                except:
                    continue
            try:
                lastAtom
            except:
                for Cindex in range(1,30):
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
            trjconvCOMMAND = '/home/osollila/Programs/gromacs/gromacs402/bin/trjconv'
            makendxCOMMAND = '/home/osollila/Programs/gromacs/gromacs402/bin/make_ndx'
        else:
            trjconvCOMMAND = 'gmx trjconv'
            makendxCOMMAND = 'gmx make_ndx'
        
        os.system('rm foo.ndx')
        os.system('echo "a ' + lastAtom + '\nq" | ' + makendxCOMMAND +' -f ' + tpr_name + ' -o foo.ndx')
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
                os.system('echo System |  ' +  trjconvCOMMAND + ' -f ' + trj_name + ' -s ' + tpr_name + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime))

            if (not os.path.isfile(xtcfoo)):
                os.system('echo "centralAtom\nSystem" |  '+  trjconvCOMMAND + ' -center -pbc mol -n foo.ndx -f ' + xtcwhole  + ' -s ' + tpr_name + ' -o ' + xtcfoo) 

            os.system('rm foo.ndx')
            os.system('rm ' + xtcwhole)
                
        #Center around the center of mass of all the g_3 carbons
        #if (not os.path.isfile(xtccentered)):        
            os.system('echo "a ' + G3atom + '\nq" | ' + makendxCOMMAND +' -f ' + tpr_name + ' -o foo.ndx')
            os.system('echo "' + G3atom + '\nSystem" |  '+  trjconvCOMMAND + ' -center -pbc mol -n foo.ndx -f ' + xtcfoo + ' -s ' + tpr_name + ' -o ' + xtccentered)
            os.system('rm ' + xtcfoo)

            
    else:
        print('Centering for other than Gromacs may not work if there are jumps over periodic boundary conditions in z-direction.')

            
            
    if download_failed == False:
        with open("system_already_run_17_11_21.out","a") as f:
            f.write(system_path+" download successfull \n")   
        if (not os.path.isfile(system_path + "/FormFactor.json")):
            try:
                if 'gromacs' in system['SOFTWARE']:
                    form_factor.FormFactor(system_path, tpr_name, xtccentered, 200, output_name,  system)
                if 'openMM' in system['SOFTWARE'] or 'NAMD' in system['SOFTWARE']:
                    form_factor.FormFactor(system_path, pdb_name, trj_name, 200, output_name,  system)
            except ValueError as e:
                #print(e)
                if "Cannot load file containing pickled data when allow_pickle=False" not in str(e):
                    raise
                #else:
                    
                #print(e)
            #    print(system_path,' Failed')
