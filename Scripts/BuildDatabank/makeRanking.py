### This code creates different types of ranking lists of simulations based on their quality against experiments.
### The ranking lists are stored in 'Data/Ranking/' folder in json format.
### The lists can be shown with the plotQuality.ipynb


import sys

sys.path.insert(1, '../BuildDatabank/')
from databankLibrary import *

databankpath = '../../Data/Simulations/'
db_data = databank(databankpath)
systems = db_data.get_systems()


#### Making list of qualities
qualities = []
for system in systems:
    quality_dict = {}
    path = databankpath + system['path']
    READMEfilepath = path + '/README.yaml'
    
    with open(READMEfilepath) as yaml_file:
        readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
    yaml_file.close()
    
    TotalQualityFilePath = path + '/SYSTEM_quality.json'
    
    FragmentQ = {}
    for lipid in readme['COMPOSITION']:
        QualityFile = path + lipid + '_FragmentQuality.json'
        #print(QualityFile)
        try:
            with open(QualityFile) as json_file:
                FragmentQ[lipid] = json.load(json_file)
            json_file.close()
        except:
            continue
        #print(lipid, FragmentQ[lipid])
    
    if (os.path.isfile(TotalQualityFilePath)):
        with open(TotalQualityFilePath) as json_file:
            FragmentQ['TotalQuality'] = json.load(json_file)
        json_file.close()
        #print(FragmentQ['TotalQuality'], '\n')

    FFQualityFilePath = path + '/FormFactorQuality.json'
    if (os.path.isfile(FFQualityFilePath)):
        with open(FFQualityFilePath) as json_file:
            FFq = json.load(json_file)
        try:
            FragmentQ['TotalQuality']
        except:
            FragmentQ['TotalQuality'] = {}
        #print(FragmentQ['TotalQuality']['headgroup'], FFq)
        try:
            FragmentQ['TotalQuality']['FFQuality'] = FFq[0]
        except:
            FragmentQ['TotalQuality']['FFQuality'] = FFq    
        json_file.close()
        #print(FragmentQ['TotalQuality'], '\n')
    
    FragmentQ['system'] = system
    qualities.append(FragmentQ)




##### Sort based on total quality of a simulation
Fragments = ['total','tails','headgroup']

for SortBasedOn in Fragments:
    #print(SortBasedOn)
    NewQualities = []
    for i in qualities:
        try:
            if i['TotalQuality'][SortBasedOn] >0:
                NewQualities.append(i)
            #totalQ = i['TotalQuality']['headgroup']
        except:
            continue

    
    SortedQualities = sorted(NewQualities, key = lambda i: i['TotalQuality'][SortBasedOn], reverse = True)
    #ShowTable(SortedQualities,'TotalQuality')
    
    outputfile = '../../Data/Ranking/SYSTEM_' + SortBasedOn + '_Ranking.json'
    with open(outputfile, "w") as fp:
        json.dump(SortedQualities, fp, default=str)
    print('Sorted based on ', SortBasedOn, ' quality and saved to', outputfile)
        

NewQualities = []
for i in qualities:
    try:
        if i['TotalQuality']['FFQuality'] >0:
            NewQualities.append(i)
            #totalQ = i['TotalQuality']['headgroup']
    except:
        continue

        
SortedQualities = sorted(NewQualities, key = lambda i: i['TotalQuality']['FFQuality'])
#ShowTable(SortedQualities,'TotalQuality')
    
outputfile = '../../Data/Ranking/SYSTEM_FormFactor_Ranking.json'
with open(outputfile, "w") as fp:
    json.dump(SortedQualities, fp, default=str)
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
                #totalQ = i['TotalQuality']['headgroup']
            except:
                continue

    
        SortedQualities = sorted(NewQualities, key = lambda i: i[lipid][SortBasedOn], reverse = True)
        #if SortedQualities:
        #    print('Quality of',SortBasedOn,' of ',lipid) #,'in simulation: ',  simulation[lipid][SortBasedOn])
        #    ShowTable(SortedQualities, lipid)
    
        if SortedQualities:
            outputfile = '../../Data/Ranking/' + lipid + '_' + SortBasedOn + '_Ranking.json'
            with open(outputfile, "w") as fp:
                json.dump(SortedQualities, fp, default=str)
            print('Quality of', SortBasedOn, ' of ' ,lipid, 'sorted and saved to', outputfile) #,'in simulation: ',  simulation[lipid][SortBasedOn])
