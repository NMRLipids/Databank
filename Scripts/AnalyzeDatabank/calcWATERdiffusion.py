#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np

import yaml
import json
import matplotlib.pyplot as plt
import mdtraj
import urllib.request
import seaborn as sns

import MDAnalysis
from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD

#from OrderParameter import *

# Download link
def download_link(doi, file):
    if "zenodo" in doi.lower():
        zenodo_entry_number = doi.split(".")[2]
        return 'https://zenodo.org/record/' + zenodo_entry_number + '/files/' + file
    else:
        print ("DOI provided: {0}".format(doi))
        print ("Repository not validated. Please upload the data for example to zenodo.org")
        return ""
    
# read mapping file
def read_mapping_file(mapping_file, atom1):
    with open(mapping_file, 'rt') as mapping_file:
            for line in mapping_file:
                if atom1 in line:
                    m_atom1 = line.split()[1]
    return m_atom1

def read_mapping_file_res(mapping_file, atom1):
    with open(mapping_file, 'rt') as mapping_file:
            for line in mapping_file:
                if atom1 in line:
                    m_res = line.split()[2]
    return m_res

def read_mapping_filePAIR(mapping_file, atom1, atom2, molname):
    with open(mapping_file, 'rt') as mapping_file:
            print(mapping_file)
            for line in mapping_file:
                if atom1 in line:
                    m_atom1 = line.split()[1]
                    try:
                        res = line.split()[2]
                    except:
                        res = molname
#                    print(m_atom1)
                if atom2 in line: 
                    m_atom2 = line.split()[1]
#                    print(m_atom2)
    return m_atom1, m_atom2, res

def make_positive_angles(x):
    for i in range(len(x)):
        if x[i] < 0:
            x[i] = np.degrees(x[i]) + 360
        else:
            x[i] = np.degrees(x[i])
    return x


# dictionary of lipids

lipids_dict = {
            'POPC' : {"REQUIRED": False,
                             "TYPE": "string",
                         },
            'POPG' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'POPS' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'POPE' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'DMPC' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'POPI' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'SAPI' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'SLPI' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'CHOL' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'DHMDMAB' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
                }


# Dictionary of other than lipid molecules.

molecules_dict = {
                
            'POT' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'SOD' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'CLA' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'CAL' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'SOL' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
                }



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
                #for molname in molecules:
                doi = readme.get('DOI')
                trj = readme.get('TRJ')
                tpr = readme.get('TPR')
                trj_name = subdir + '/' + readme.get('TRJ')[0][0]
                tpr_name = subdir + '/' + readme.get('TPR')[0][0]
                gro_name = subdir + '/conf.gro'
                trj_url = download_link(doi, trj[0][0])
                tpr_url = download_link(doi, tpr[0][0])
                EQtime=float(readme.get('TIMELEFTOUT'))*1000

                densityFOLDERS = subdir.replace("Simulations","WATERdiffusion")    
                outfilename = str(densityFOLDERS) + '/WATERlateralMSD.xvg'
                #indexfilename = str(densityFOLDERS) + '/' + molname + 'index.ndx' 

                #    try:
                #        nmol = int(readme['N' + molname])
                #    except:
                #        try:
                #            nmol = int(readme['N' + molname][0])
                #        except:
                #            continue

                print(outfilename)                    
                if not os.path.isfile(outfilename):
                                            
                    print('Analyzing in '+filepath)
                        
                    #Download tpr and xtc files to same directory where dictionary and data are located
                    if (not os.path.isfile(tpr_name)):
                        response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
                    if (not os.path.isfile(trj_name)):
                        response = urllib.request.urlretrieve(trj_url, trj_name)
                        
                    fig= plt.figure(figsize=(12,9))
                    if (not os.path.isfile(gro_name)):
                        get_ipython().system('echo System | gmx trjconv -f {trj_name} -s {tpr_name}  -dump 0 -o {gro_name}')
                        
                    xtcwhole=subdir + '/whole.xtc'
                    if (not os.path.isfile(xtcwhole)):
                        get_ipython().system('echo System | gmx trjconv -f {trj_name} -s {tpr_name} -o {xtcwhole} -pbc mol -b {EQtime}')
                        
                        
                    os.system('mkdir -p ' + densityFOLDERS)
                    os.system('cp ' + READMEfilepath + ' ' + densityFOLDERS)

                    #u = MDAnalysis.Universe(gro_name, trj_name)
                    water_mapping_file = '../BuildDatabank/mapping_files/'+readme['MAPPING_DICT']['SOL']
                    waterO = read_mapping_file(water_mapping_file,'M_O_M')
                    #waterRES = readme.get('SOL')
                    
                    #os.system('echo System | gmx trjconv -f ' + xtc + ' -s ' + tpr + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime))

                    ndxSTR = 'keep 0 \na ' + waterO  + ' \nkeep 1 \nq'
                    #print(ndxSTR)
                    os.system('echo "' + ndxSTR + '" | gmx make_ndx -f ' + tpr_name + ' -o ' + densityFOLDERS + '/SOLindex.ndx') 
                    os.system('gmx msd -f ' + xtcwhole + ' -s ' + tpr_name + ' -n  ' + densityFOLDERS + '/SOLindex.ndx -o ' + outfilename + ' -trestart 1000 -lateral z')
                    
                    #select = "byres name " + waterOX
                    #MSD_analysis = MSD(u, select, 0, 1000, 20)
                    #MSD_analysis.run()
                
                    #NcenterMOL = 0
                    #for lipid in lipids_dict:
                    #    if readme.get('N' + lipid):
                    #        if int(readme.get('N' + lipid)[0]) > NcenterMOL:
                    #            mapping_file = '../BuildDatabank/mapping_files/'+readme['MAPPING_DICT'][lipid] # readme.get('MAPPING')[0][0]
                    #            try:
                    #                centerMOL = read_mapping_file_res(mapping_file,'M_G3P2_M')
                    #            except:
                    #                centerMOL = readme.get(lipid)
                    #            NcenterMOL = int(readme.get('N' + lipid)[0])

                    
                        #if (not os.path.isfile(outfilename)):
                        #    try:
                        #        mapping_file = '../BuildDatabank/mapping_files/'+readme['MAPPING_DICT'][molname] # readme.get('MAPPING')[0][0]
                        #        resname = read_mapping_file_res(mapping_file,'M_G3P2_M')
                        #    except:
                        #        resname = readme.get(molname)
                        #    try:
                        #        get_ipython().system('echo "r {centerMOL} \n name 0 centerMOL \n r {resname} \n q" | gmx make_ndx -f {tpr_name} -o {indexfilename} -n /media/osollila/Data/NMRlipids/NMRlipidsIVPEandPG/scripts/init.ndx')
                        #        get_ipython().system('echo "centerMOL \n  {resname} \n q" | gmx density -f {xtcwhole} -s {tpr_name} -o {outfilename} -b {EQtime} -n {indexfilename} -xvg none -center -dens number')
                        #    except:
                        #        print('Density calculation failed.')

            







