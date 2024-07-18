#!/usr/bin/env python3
# coding: utf-8

import sys, os
import numpy as np
import json
from tqdm import tqdm
import traceback

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from DatabankLib.databankLibrary import initialize_databank
from DatabankLib import NMLDB_SIMU_PATH

systems = initialize_databank()

cErr = 0
cUpp = 0

for system in tqdm(systems, desc='Scan over Databank systems'):
    thickFN = os.path.join(NMLDB_SIMU_PATH, system['path'], 'thickness.json')
    if (os.path.isfile(thickFN)):
        # file exist. Skipping
        continue

    WaterDensity_name = os.path.join(NMLDB_SIMU_PATH, system['path'], 'WaterDensity.json')
    LipidDensity_name = os.path.join(NMLDB_SIMU_PATH, system['path'], 'LipidDensity.json')
    print(LipidDensity_name)
    try:
        with open(WaterDensity_name) as f:
            WaterDensity = json.load(f)
        with open(LipidDensity_name) as f:
            LipidDensity = json.load(f)
        wd = np.array(WaterDensity)
        ld = np.array(LipidDensity)
        idx = np.argwhere(np.diff(np.sign(wd[:,1] - ld[:,1]))).flatten()
        if len(idx) < 2:
            print("Dehydrated sample! Standard thickness rule doesn't work. Will extract boxZ.")
            thickness = wd[-1, 0] - wd[0, 0]
        else:
            thickness = wd[idx[1], 0] - wd[idx[0], 0]
        with open(thickFN, 'w') as f:
            json.dump(thickness, f)
            cUpp += 1
    except Exception as e:
        print('Calculation failed for ' +  system['path'])
        print(str(e))
        print(traceback.format_exc())
        cErr += 1

print(f"""
===========================================
Summary:
 Thickness updated for: {cUpp} systems
 Thickness calculation failed: {cErr} systems
""")