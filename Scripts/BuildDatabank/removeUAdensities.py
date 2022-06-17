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
    try:
        if system['UNITEDATOM_DICT']:
            print(system['path'])
            os.system('rm ' + system['path'] + '/*Density*')
            os.system('rm ' + system['path'] + '/FormFactor*')
            #os.system('rm ' + system['path'] + '/FormFactorQuality*')
            #os.system('ls ' + system['path'] + '/thickness.json')
    except:
        continue
