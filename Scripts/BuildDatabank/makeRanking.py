#!/usr/bin/env python3
# coding: utf-8
"""
`makeRanking.py` creates different types of ranking lists of simulations based on their quality against experiments.
The ranking lists are stored in `Data/Ranking/` folder in json format.
The lists can be shown with the `plotQuality.ipynb`
"""

import sys, os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from DatabankLib.databankLibrary import *
from DatabankLib.jsonEncoders import CompactJSONEncoder

systems = initialize_databank()

#### Making list of qualities
qualities = []
for system in systems:
    quality_dict = {}
    path = os.path.join(NMLDB_SIMU_PATH, system['path'])
    READMEfilepath = os.path.join(path, 'README.yaml')
    
    with open(READMEfilepath, 'r') as yaml_file:
        readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
    
    TotalQualityFilePath = os.path.join(path, 'SYSTEM_quality.json')
    
    FragmentQ = {}
    for lipid in readme['COMPOSITION']:
        QualityFile = os.path.join(path, lipid + '_FragmentQuality.json')
        try:
            with open(QualityFile, 'r') as json_file:
                FragmentQ[lipid] = json.load(json_file)
        except:
            continue
    
    if (os.path.isfile(TotalQualityFilePath)):
        with open(TotalQualityFilePath, 'r') as json_file:
            FragmentQ['TotalQuality'] = json.load(json_file)

    FFQualityFilePath = os.path.join(path, 'FormFactorQuality.json')
    if (os.path.isfile(FFQualityFilePath)):
        with open(FFQualityFilePath) as json_file:
            FFq = json.load(json_file)
        try:
            FragmentQ['TotalQuality']
        except:
            FragmentQ['TotalQuality'] = {}
        try:
            FragmentQ['TotalQuality']['FFQuality'] = FFq[0]
        except:
            FragmentQ['TotalQuality']['FFQuality'] = FFq    
        json_file.close()
    
    FragmentQ['system'] = system
    qualities.append(FragmentQ)


##### Sort based on total quality of a simulation
Fragments = ['total','tails','headgroup']

for SortBasedOn in Fragments:
    NewQualities = []
    for i in qualities:
        try:
            if i['TotalQuality'][SortBasedOn] >0:
                NewQualities.append(i)
        except:
            continue
    
    SortedQualities = sorted(NewQualities, key = lambda i: i['TotalQuality'][SortBasedOn], reverse = True)
    
    outputfile = os.path.join(NMLDB_DATA_PATH, 'Ranking', 'SYSTEM_' + SortBasedOn + '_Ranking.json')
    with open(outputfile, "w") as fp:
        json.dump(SortedQualities, fp, default=str, cls=CompactJSONEncoder)
    print('Sorted based on ', SortBasedOn, ' quality and saved to', outputfile)
        

NewQualities = []
for i in qualities:
    try:
        if i['TotalQuality']['FFQuality'] >0:
            NewQualities.append(i)
    except:
        continue

        
SortedQualities = sorted(NewQualities, key = lambda i: i['TotalQuality']['FFQuality'])
    
outputfile = os.path.join(NMLDB_DATA_PATH, 'Ranking', 'SYSTEM_FormFactor_Ranking.json')
with open(outputfile, "w") as fp:
    json.dump(SortedQualities, fp, default=str, cls=CompactJSONEncoder)
print('Sorted based on form factor quality and saved to', outputfile)    


##### Sorting best simulations for each lipid
Fragments = ['total','sn-1','sn-2','headgroup']

for SortBasedOn in Fragments:
    for lipid in lipids_dict:
        NewQualities = []
        for i in qualities:
            try:
                if i[lipid][SortBasedOn] >0:
                    NewQualities.append(i)
            except:
                continue
    
        SortedQualities = sorted(NewQualities, key = lambda i: i[lipid][SortBasedOn], reverse = True)
    
        if SortedQualities:
            outputfile = os.path.join(NMLDB_DATA_PATH, 'Ranking', lipid + '_' + SortBasedOn + '_Ranking.json')
            with open(outputfile, "w") as fp:
                json.dump(SortedQualities, fp, default=str, cls=CompactJSONEncoder)
            print('Quality of', SortBasedOn, ' of ' ,lipid, 'sorted and saved to', outputfile) 
            #,'in simulation: ',  simulation[lipid][SortBasedOn])
