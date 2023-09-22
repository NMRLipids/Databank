import os
import sys
import numpy as np
import json
import matplotlib.pyplot as plt
import MDAnalysis
import urllib.request
import yaml

#sys.path.insert(1, '../BuildDatabank/')
#from databankLibrary import download_link, lipids_dict, databank

#path = '../../Data/Simulations/'
#db_data = databank(path)
#systems = db_data.get_systems()


databankPath =  '../../'
sys.path.insert(1, databankPath + '/Scripts/BuildDatabank/')
from databankLibrary import * 
systems = initialize_databank(databankPath)


IDs = []
for system in systems:
    if 'ID' in system.keys():
        IDs.append(system['ID'])

print('Adding IDs to the following to these systems:')
for system in systems:
    if 'ID' not in system.keys():
        NewID = max(IDs) + 1
        system['ID'] = NewID
        READMEpath = '../../Data/Simulations/' + system['path'] + 'README.yaml'
        #print(READMEpath)
        system.pop('path',None)
        with open(READMEpath, 'w') as f:
            yaml.dump(system,f, sort_keys=False)
        IDs.append(system['ID'])
        print(READMEpath, NewID)
        
print('Largest ID: ', max(IDs))

#ind = 1
#for system in systems:
#    system['ID'] = ind
#    ind += 1
#    READMEpath = system['path'] + 'README.yaml'
#    system.pop('path',None)
#    with open(READMEpath, 'w') as f:
#        yaml.dump(system,f, sort_keys=False)
#    #print(READMEpath)
        

