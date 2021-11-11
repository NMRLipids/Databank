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
import form_factor


path = '../../Data/Simulations/'

db_data = databank(path)
systems = db_data.get_systems()

for system in systems:
    #download trajectory and gro files
    system_path = system['path']
    
    output_name = "form_factors"
    
    trj_name = system['path'] + system['TRJ'][0][0]
    tpr_name = system['path'] + system['TPR'][0][0]
    trj_url = download_link(system['DOI'], system['TRJ'][0][0])
    tpr_url = download_link(system['DOI'], system['TPR'][0][0])

    if (not os.path.isfile(tpr_name)):
        response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
    if (not os.path.isfile(trj_name)):
        response = urllib.request.urlretrieve(trj_url, trj_name)
        
    form_factor.FormFactor(system_path, tpr_name, trj_name, 200, output_name, 'resname POPS', system)
