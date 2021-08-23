#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import numpy as np

import yaml
import json
import matplotlib.pyplot as plt
import mdtraj
import urllib.request
import seaborn as sns

sys.path.insert(1, '../BuildDatabank/')
from databankLibrary import download_link, read_mapping_file, read_mapping_file_res, read_mapping_filePAIR, make_positive_angles, lipids_dict, molecules_dict

# Calculate density profiles of all molecules

#molecules = {'POPS','POPE','POPG','POPC','POT','SOD','CLA','CAL','SOL'}
molecules = []
for key in lipids_dict:
    molecules.append(key)
for key in molecules_dict:
    molecules.append(key)
#lipids = {'POPS','POPE','POPG','POPC'}
#lipids = {'POPE'}
#atom1 = 'M_G3P2_M'
#atom2 = 'M_G3C5O1_M'   #for POPG
#atom2 = 'M_G3N6_M'

colors = {'POPC' :'black','POPS':'red','POPE':'blue','POPG':'green'}

h = []


for subdir, dirs, files in os.walk(r'../../Data/Simulations/'):
    for filename in files:
        filepath = subdir + os.sep + filename
        if filepath.endswith("README.yaml"):
            READMEfilepath = subdir + '/README.yaml'
            with open(READMEfilepath) as yaml_file:
                readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
                molname = 'SOL' 
                doi = readme.get('DOI')
                trj = readme.get('TRJ')
                tpr = readme.get('TPR')
                trj_name = subdir + '/' + readme.get('TRJ')[0][0]
                tpr_name = subdir + '/' + readme.get('TPR')[0][0]
                gro_name = subdir + '/conf.gro'
                trj_url = download_link(doi, trj[0][0])
                tpr_url = download_link(doi, tpr[0][0])
                EQtime=float(readme.get('TIMELEFTOUT'))*1000
                
                outFOLDERS = subdir.replace("Simulations","WATERT1")    
                CorrFname = str(outFOLDERS) + '/' + molname + 'waterCORRF.xvg'
                indexfilename = str(outFOLDERS) + '/' + molname + 'index.ndx' 

                print(molname, readme.get('N' + molname), CorrFname)                    
                if not os.path.isfile(CorrFname):
                    print('Analyzing '+molname+' in '+filepath)
                        
                    #Download tpr and xtc files to same directory where dictionary and data are located
                    if (not os.path.isfile(tpr_name)):
                        response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
                    if (not os.path.isfile(trj_name)):
                        response = urllib.request.urlretrieve(trj_url, trj_name)
                        
                    #fig= plt.figure(figsize=(12,9))
                    if (not os.path.isfile(gro_name)):
                        os.system('echo System | gmx trjconv -f ' + trj_name + ' -s ' + tpr_name + ' -dump 0 -o '  + gro_name)
                        
                    xtcwhole=subdir + '/whole.xtc'
                    if (not os.path.isfile(xtcwhole)):
                        os.system('echo System | gmx trjconv -f ' + trj_name + ' -s ' + tpr_name + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime))
                        
                        
                    os.system('mkdir -p ' + outFOLDERS)
                    os.system('cp ' + READMEfilepath + ' ' + outFOLDERS)

                    mapping_file = '../BuildDatabank/mapping_files/'+readme['MAPPING_DICT'][molname]
                    resname = readme.get(molname)
                    Oatom = read_mapping_file(mapping_file, 'M_O_M')
                    H1atom = read_mapping_file(mapping_file, 'M_H1_M')
                    H2atom = read_mapping_file(mapping_file, 'M_H2_M')
                    #try:
                    os.system('echo "r ' + resname + ' & a ' + H1atom + ' | a ' + Oatom + ' \n q" | gmx make_ndx -f ' + tpr_name + ' -o ' + indexfilename + ' -n /media/osollila/Data/NMRlipids/NMRlipidsIVPEandPG/scripts/init.ndx')
                    os.system('gmx rotacf -f ' + xtcwhole + ' -s ' + tpr_name + ' -o ' + CorrFname + ' -b ' + str(EQtime) + ' -n ' + indexfilename + ' -xvg none -P 2 -d')
                    #except:
                    #    print('Calculation failed.')

            







