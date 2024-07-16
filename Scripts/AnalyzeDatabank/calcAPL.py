#!/usr/bin/env python3
# coding: utf-8

import os, sys
import json

sys.path.append('..')
from DatabankLib import NMLDB_SIMU_PATH
from DatabankLib import jsonEncoders
from DatabankLib.databankLibrary import *

### Initializing the databank
systems = initialize_databank()

### Loop over simulations in the databank
for system in systems:
    ## reading software and file path for simulations
    software = system['SOFTWARE']
    path = system['path']

    ## this is checking if area per lipid is already calculated for the systems
    outfilename = os.path.join(NMLDB_SIMU_PATH,  path, 'apl.json')
    if os.path.isfile(outfilename):
        continue
    
    print('Analyzing: ', path)
    print('Will write into: ', outfilename)

    ## calculates the total number of lipids
    Nlipid = GetNlipids(system)

    ## makes MDAnalysis universe from the system. This also downloads the data if not yet locally available
    u = system2MDanalysisUniverse(system)

    if u is None:
        print('Generation of MDAnalysis universe failed in folder', path)
        continue
    
    ## this calculates the area per lipid as a function of time and stores it in the databank
    apl = {}
    for ts in tqdm(u.trajectory, desc='Scanning the trajectory'):
        if u.trajectory.time >= system['TIMELEFTOUT']*1000:
            dimensions = u.dimensions
            aplFrame = u.dimensions[0]*u.dimensions[1]*2/Nlipid
            apl[u.trajectory.time] = aplFrame

    with open(outfilename, 'w') as f:
        json.dump(apl, f, cls=jsonEncoders.CompactJSONEncoder)
        
