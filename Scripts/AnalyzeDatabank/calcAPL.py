import os
import sys
import numpy as np
import json
import matplotlib.pyplot as plt
import MDAnalysis
import urllib.request
import yaml

sys.path.insert(1, '../BuildDatabank/')
from databankLibrary import download_link, lipids_dict, molecule_numbers_dict,databank

path = '../../Data/Simulations/'
db_data = databank(path)
systems = db_data.get_systems()

for system in systems:
    Nlipid = 0
    path = system['path']
    outfilename = path + 'apl.json'
    if os.path.isfile(outfilename):
        continue

    print('Analyzing: ', system['path'])

    for molecule in system['COMPOSITION']:
        if molecule in lipids_dict:
            Nlipid += np.sum(system['COMPOSITION'][molecule]['COUNT'])

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
                        
    if (not os.path.isfile(trj_name)):
        response = urllib.request.urlretrieve(trj_url, trj_name)

    try:
        u = MDAnalysis.Universe(tpr_name, trj_name)
    except:
        conf = path + '/conf.gro'
        print("Generating conf.gro because MDAnalysis cannot read tpr version")
        os.system('echo System | gmx trjconv -s '+ tpr_name + ' -f '+ trj_name + ' -dump 0 -o ' + conf)
        u = MDAnalysis.Universe(conf, trj_name)

    apl = {}
    for ts in u.trajectory:
        if u.trajectory.time >= system['TIMELEFTOUT']*1000:
            dimensions = u.dimensions
            #print(dimensions)
            aplFrame = u.dimensions[0]*u.dimensions[1]*2/Nlipid
            apl[u.trajectory.time] = aplFrame
            #print(apl)

    with open(outfilename, 'w') as f:
        json.dump(apl,f)
        
    #print(outfilename)
        

