#!/usr/bin/env python3
# coding: utf-8

import sys, os
import numpy as np
import json
from tqdm import tqdm

sys.path.append('..')
from DatabankLib.databankLibrary import databank
from DatabankLib import NMLDB_SIMU_PATH

db_data = databank()
systems = db_data.get_systems()

for system in tqdm(systems, desc='Scan over Databank systems'):
    thickFN = os.path.join(NMLDB_SIMU_PATH, system['path'], 'thickness.json')
    if (os.path.isfile(thickFN)):
        # file exist. Skipping
        continue

    WaterDensity_name = os.path.join(NMLDB_SIMU_PATH, system['path'], 'WaterDensity.json')
    LipidDensity_name = os.path.join(NMLDB_SIMU_PATH, system['path'], 'LipidDensity.json')
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
        f = open(thickFN, 'w')
        json.dump(thickness,f)
    except:
        print('Calculation failed for ' +  system['path'])
