import os
import sys
import numpy as np
import json
import matplotlib.pyplot as plt
import MDAnalysis
import urllib.request
import yaml

sys.path.insert(1, '../BuildDatabank/')
from databankLibrary import download_link, lipids_dict, databank

path = '../../Data/Simulations/'
db_data = databank(path)
systems = db_data.get_systems()

for system in systems:
    Nlipid = 0
    path = system['path']
    outfilename = path + 'conf.pdb'
    if os.path.isfile(outfilename):
        continue

    print('Analyzing: ', system['path'])

    doi = system.get('DOI')
    trj = system.get('TRJ')
    tpr = system.get('TPR')
    trj_name = path + system.get('TRJ')[0][0]
    tpr_name = path + system.get('TPR')[0][0]
    trj_url = download_link(doi, trj[0][0])
    tpr_url = download_link(doi, tpr[0][0])
    #print(trj_name,tpr_name)
    
    if (not os.path.isfile(tpr_name)):
        response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
    #if (not os.path.isfile(trj_name)):
    #    response = urllib.request.urlretrieve(trj_url, trj_name)

    if 'gromacs' in system['SOFTWARE']:
        os.system('gmx editconf -f '+ tpr_name + ' -o ' + outfilename)


        

