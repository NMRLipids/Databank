#Lipid flip-flops
#loops over trajectories and checks if some molecule (lipid or cholesterol) in membrane flips
import os
import sys
import yaml
import json
import matplotlib.pyplot as plt
import numpy as np
import math
from random import randint

from matplotlib import cm
from scipy.stats import norm

import urllib.request
from urllib.error import URLError,HTTPError,ContentTooShortError

import MDAnalysis as mda

sys.path.insert(1, '../BuildDatabank/')
from databankLibrary import download_link, lipids_dict

def getLipids(readme):
        lipids = []
        for key in lipids_dict.keys():
            try:
                if readme['N'+key] != [0,0]: 
                    lipids.append(key)
            except KeyError:
                continue
        return lipids
        
def getHeadgroup(readme, lipid):
    headgroup = ""
    m_file = readme['MAPPING_DICT'][lipid]
    with open('../BuildDatabank/mapping_files/'+m_file,"r") as f:
        for line in f:
            if 'M_G3' in line:
                headgroup = headgroup + ' ' + line.split()[1]
  
    return headgroup

def membraneCentreOfMass(universe, readme):
    lipids = []
    for key_mol in lipids_dict:
        selection = ""
        if key_mol in readme['MAPPING_DICT'].keys():
            m_file = readme['MAPPING_DICT'][key_mol]
            with open('../BuildDatabank/mapping_files/'+m_file,"r") as f:
                for line in f:
                    if len(line.split()) > 2 and "Individual atoms" not in line:
                        selection = selection + "(resname " + line.split()[2] + " and name " + line.split()[1] + ") or "
                    elif "Individual atoms" in line:
                       continue
                    else:
                        selection = "resname " + readme[key_mol]
                        #print(selection)
                        break
        selection = selection.rstrip(' or ')
        #print("selection    " + selection)
        molecules = universe.select_atoms(selection)
        #print("molecules")
        #print(molecules)
        if molecules.n_residues > 0:
            lipids.append(universe.select_atoms(selection))
            #print(lipids) 
    # join all the selected the lipids together to make a selection of the entire membrane and calculate the
    # z component of the centre of mass of the membrane
    membrane = universe.select_atoms("")
    R_membrane_z = 0
    if lipids!= []:
        for i in range(0,len(lipids)):
            membrane = membrane + lipids[i]
            #print("membrane") 
            #print(membrane)  
    R_membrane = membrane.center_of_mass()[2]
    return R_membrane
    
# check if lipid headgroup switches to other side of the centre of mass of the membrane 


    
    


#if (not os.path.isdir('../../Data/Simulations/')): 
#    os.mkdir('../../Data/Simulations/')

for subdir, dirs, files in os.walk(r'../../Data/Simulations/'): 
    for filename1 in files:
        filepath = subdir + os.sep + filename1        
        if filepath.endswith("README.yaml"):
            READMEfilepathSimulation = subdir + '/README.yaml'
            with open(READMEfilepathSimulation) as yaml_file_sim:
                #print('toimii')
                readme = yaml.load(yaml_file_sim, Loader=yaml.FullLoader)
                indexingPath = "/".join(filepath.split("/")[4:8])
                print(indexingPath)
                DATAdir = '../../../../tst/FlipFlop/' + indexingPath + '/'                                     #change directory when code works properly!!!
                
                doi = readme['DOI']
                trj=readme['TRJ'][0][0]
                tpr=readme['TPR'][0][0]
                trj_name = subdir + '/' + trj
                tpr_name = subdir + '/' + tpr
                trj_url = download_link(doi, trj)
                tpr_url = download_link(doi, tpr)
                EQtime=float(readme['TIMELEFTOUT'])*1000
                
                #Download tpr and xtc files to a tmp directory in dir wrk defined in readme file / or should these be put into the same directory as README.yaml ???
                if (not os.path.isfile(tpr_name)):
                    response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
                if (not os.path.isfile(trj_name)):
                    response = urllib.request.urlretrieve(trj_url, trj_name)
                    
                xtcwhole=subdir + '/whole.xtc'
                if (not os.path.isfile(xtcwhole)):
                    os.system('echo System | gmx trjconv -f ' + trj_name + ' -s ' + tpr_name + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime))
                    
                u = mda.Universe(tpr_name,xtcwhole)
                
                lipids = getLipids(readme)
                
                end = len(u.trajectory)
                
                dt = u.trajectory.dt
                start = int(EQtime*dt)
                flipflops = 0
                time_diff = 500
                skip = int(time_diff / dt)
                
                frames_lipids = []
                
                for lipid in lipids:
                    previous_leaflet = []
                    #every 500th frame for checking flip flops
                    for ts in u.trajectory[start:end:skip]:
                        print('frame ' + str(ts.frame))
                        headgroup = getHeadgroup(readme, lipid)
                       # print(headgroup)

                        hg_selection = u.select_atoms('resname ' + lipid + ' and name' + headgroup).split('residue')
                        # ei toimi
                        R_m = membraneCentreOfMass(u, readme)
                        hg_leaflet_fr = []
                        
                        for hg in hg_selection:  
                            hg_z = hg.center_of_mass()[2]
                            if hg_z < R_m:
                                hg_leaflet_fr.append('l1') #in leaflet 1
                            if hg_z > R_m:
                                hg_leaflet_fr.append('l2') #in leaflet 2
                               
                        
                        #check if lipid has changed leaflet between frames
                        if ts.frame != 0:
                            for i in range(0,len(hg_leaflet_fr)):
                                if hg_leaflet_fr[i] != previous_leaflet[i]:
                                    print("FLIPFLOP !!!!!!")
                                    flipflops += 1
                                    frames_lipids.append([ts.frame + '   ' + lipid])
                                    
                                    # save frame nr and lipid name
                                    with open(flipflop_file, 'w') as f:
                                        for fr in frames_lipids:
                                            f.write(fr + '\n')
                                    f.close()
                                    
                        
                        previous_leaflet = hg_leaflet_fr.copy()    
                            
                        
                       
                            
                if flipflops > 0:
                    print('number of flip flops in trajectory is ' + str(flipflops) + " for simulation in directory " + indexingPath)
                    flipflop_file = '../../../../tst/FlipFlop/' + indexingPath + '/flipflop.txt'
                    with open(flipflop_file, 'w') as f:
                        f.write('# frame   lipid\n')
                        
                    f.close()
                else:
                    print("no flip flops in trajectory")
                    
                
                
