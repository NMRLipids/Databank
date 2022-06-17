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
    WaterDensity_name = system['path'] + 'WaterDensity.json'
    try:
        f = open(WaterDensity_name)
        #print('Density file not found')
        WaterDensity = json.load(f)
        center = round(len(WaterDensity)/2)
        #print(WaterDensity[center][1], WaterDensity[0][1])
        if WaterDensity[center][1] > WaterDensity[0][1]:
            print(system['path'])
            system['WARNINGS'] = {}
            system['WARNINGS']['PBC'] = 'z-jumps'
            #print(system)

            outfileDICT= system['path'] + '/README.yaml'
            with open(outfileDICT, 'w') as f:
                yaml.dump(system,f, sort_keys=False)
            
            #os.system('rm ' + system['path'] + '/*Density*')
            #os.system('rm ' + system['path'] + '/FormFactor*')
    except:
        pass
        #print('Density file not found from ' +  system['path'])
