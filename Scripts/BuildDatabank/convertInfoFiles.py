#convert README.yaml to new format where COMPOSITION contains molecule names, mapping file names and molecule numbers

import yaml
import os

def convertFile(info_file):
    new_info_file = info_file.copy()
    del new_info_file['MAPPING_DICT']
    new_info_file['COMPOSITION'] = {}
    for key in info_file['MAPPING_DICT'].keys():
        mol_name = info_file[key]
        mapping_name = info_file['MAPPING_DICT'][key]
        new_info_file['COMPOSITION'][key] = {}
        new_info_file['COMPOSITION'][key]['NAME'] = mol_name
        new_info_file['COMPOSITION'][key]['MAPPING'] = mapping_name
        del new_info_file[key]
    new_info_file['DIR_WRK'] = '/media/osollila/Data/tmp/DATABANK/'
    return new_info_file
    
#def addWRKdir(info_file):
#    new_info_file = info_file.copy()
#    new_info_file['DIR_WRK'] = '/media/akiirikk/DATADRIVE1/tietokanta/Data/tmp/DATABANK'
#    return new_info_file


for subdir, dirs, files in os.walk(r'./info_files/'):
    for filename in files:
        filepath = subdir + os.sep + filename
        if filename.endswith('.yaml'):
            new_info_file = {}
            with open(filepath) as f:
                info_file = yaml.load(f, Loader=yaml.FullLoader)
                try:
                    tst = info_file['COMPOSITION']
                    continue
                except:
                    print('Converting', info_file['DOI'])
                new_info_file = convertFile(info_file)
                print(new_info_file)
            f.close()
            #info_file_path = subdir + os.sep + 'tst.yaml'
            with open(filepath, 'w') as f2:
                yaml.dump(new_info_file,f2, sort_keys=False)
            f2.close()
