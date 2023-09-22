import os
import sys
import numpy as np
import json
import matplotlib.pyplot as plt
import MDAnalysis
import urllib.request
import yaml

### Initializing the databank

### This is the NMRlipids databank repository path
databankPath = '../../'

sys.path.insert(1, databankPath + '/Scripts/BuildDatabank/')
from databankLibrary import *


systems = initialize_databank(databankPath)


### Loop over simulations in the databank
for system in systems:
    ## reading software and file path for simulations
    software=system['SOFTWARE']
    path = system['path']

    ## this is checking if area per lipid is already calculated for the systems
    outfilename = databankPath + '/Data/Simulations/' +  path + 'apl.json'
    if os.path.isfile(outfilename):
        continue

    print('Analyzing: ', path)

    ## calculates the total number of lipids
    Nlipid = GetNlipids(system)

    ## makes MDAnalysis universe from the system. This also downloads the data if not yet locally available
    u = system2MDanalysisUniverse(system)

    if u is None:
        print('Generation of MDAnalysis universe failed in folder', path)
        continue
    
    ## this calculates the area per lipid as a function of time and stores it in the databank
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
        

