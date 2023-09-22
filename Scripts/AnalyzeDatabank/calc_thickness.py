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

path = '../../'

db_data = databank(path)
systems = db_data.get_systems()

for system in systems:
    WaterDensity_name = '../../' + system['path'] + 'WaterDensity.json'
    LipidDensity_name = '../../' + system['path'] + 'LipidDensity.json'
    print(LipidDensity_name)
    try:
        f = open(WaterDensity_name)
        WaterDensity = json.load(f)
        f = open(LipidDensity_name)
        LipidDensity = json.load(f)
        wd=np.empty(len(WaterDensity))
        ld=np.empty(len(WaterDensity))
        for i in range(len(WaterDensity)):
            wd[i] = WaterDensity[i][1]
            ld[i] = LipidDensity[i][1]
        idx = np.argwhere(np.diff(np.sign(wd - ld ))).flatten()
        #print(system['path'])
        thickness = WaterDensity[idx[1]][0] - WaterDensity[idx[0]][0]
        f = open('../../' + system['path'] + '/thickness.json', 'w')
        json.dump(thickness,f)
    except:
        print('Calculation failed for ' +  system['path'])
