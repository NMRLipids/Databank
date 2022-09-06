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

ind = 1
for system in systems:
    system['ID'] = ind
    ind += 1
    READMEpath = system['path'] + 'README.yaml'
    system.pop('path',None)
    with open(READMEpath, 'w') as f:
        yaml.dump(system,f, sort_keys=False)
    #print(READMEpath)
        

