import os
import sys
import numpy as np
import json
import matplotlib.pyplot as plt
import MDAnalysis as mda
import urllib.request
import yaml

import urllib.request
from urllib.error import URLError,HTTPError,ContentTooShortError

# From time monitoring
from tqdm import tqdm

import socket

sys.path.insert(1, '../BuildDatabank/')
from databankLibrary import download_link, lipids_dict, databank

path = '../../Data/Simulations/'

db_data = databank(path)
systems = db_data.get_systems()

for system in systems:
    trj_name = system['path'] + system['TRJ'][0][0]
    tpr_name = system['path'] + system['TPR'][0][0]
    u = mda.Universe(tpr_name,trj_name)
    for molecule in system['COMPOSITION']:
        #print(system['COMPOSITION'][molecule]['MAPPING'])
        m_file = system['COMPOSITION'][molecule]['MAPPING']
        with open('./mapping_files/'+m_file,"r") as f:
            for line in f:
                #if len(line.split()) > 2 and "Individual atoms" not in line:
                #    selection = 'resname ' + line.split()[2] + ' and name ' +  line.split()[1] 
                #el
                selection = ''
                try:
                    if system['UNITEDATOM_DICT'] and 'H' in line.split()[0]:
                        continue
                    elif "Individual atoms" in line:
                        continue
                    elif len(line.split()) > 2 and "Individual atoms" not in line:
                        selection = 'resname ' + line.split()[2] + ' and name ' +  line.split()[1] 
                    else:
                        selection = 'resname ' + system['COMPOSITION'][molecule]['NAME'] + ' and name ' +  line.split()[1]
                except:
                    if "Individual atoms" in line:
                        continue
                    elif len(line.split()) > 2:
                        selection = 'resname ' + line.split()[2] + ' and name ' +  line.split()[1] 
                    else:
                        selection = 'resname ' + system['COMPOSITION'][molecule]['NAME'] + ' and name ' +  line.split()[1]

                NatomsFromMapping = len(u.select_atoms(selection))
                NatomsFromREADME = np.sum(system['COMPOSITION'][molecule]['COUNT'])
                #print(len(u.select_atoms(selection)), np.sum(system['COMPOSITION'][molecule]['COUNT']))
                if NatomsFromMapping !=  NatomsFromREADME and molecule not in lipids_dict:
                    print(system['path'])
                    #os.system('rm ' + system['path'] + '/*Density*')
                    #os.system('rm ' + system['path'] + '/FormFactor*')
                    print(system['COMPOSITION'][molecule]['NAME'])
                    print(selection)
                    print(NatomsFromMapping, NatomsFromREADME)
