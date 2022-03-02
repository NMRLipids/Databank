#transform dat files to json

import os
import numpy as np

import yaml
import json
import sys

data_file = sys.argv[1]
print(data_file)
if '.dat' in data_file:
    outfile = data_file.replace(".dat","_FormFactor.json")
if '.xff' in data_file:
    outfile = data_file.replace(".xff","_FormFactor.json")
print(outfile)


#data_file ="/media/akiirikk/DATADRIVE1/tietokanta/NMRLipids_Databank/Databank/Data/experiments/FormFactors/10.1016/j.bbamem.2011.07.022/17/DSPC_ULV_60Cin0D.xff"

#save_path = "/media/akiirikk/DATADRIVE1/tietokanta/NMRLipids_Databank/Databank/Data/experiments/FormFactors/10.1016/j.bbamem.2011.07.022/17/DSPC_ULV_60Cin0D"



#outfile = save_path + ".json"
data = []

with open(data_file) as OPfile:
    lines=OPfile.readlines()
    for line in lines:
        if "#" in line or line == "\n":
            continue
        print(line.split())
        if len(line.split()) == 3:
            values = [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])] #, float(line.split()[3])]
        # SAMULI: If error is not given in the original file, error of 0.02 is assumed.
        if len(line.split()) == 2:
            values = [float(line.split()[0]), float(line.split()[1]), 0.02] #, float(line.split()[3])]

            
        data.append(values)
        
        with open(outfile, 'w') as f:
            json.dump(data,f)

