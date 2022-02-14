#convert mapping files from txt file to dictionary format

import os
import sys
import yaml



for filename in os.scandir('./mapping_files/'):
    if filename.path.endswith('mappingDESCRIBE.txt') or filename.path.endswith('mappingPOPCulmschneider.txt') or filename.path == './mapping_files/mapping_files_yaml':
        continue
#    elif filename.path.endswith('mappingPOPClipid17ecc.txt'):
    #    if filename == 'mappingPOPCslipids.txt':
  #          print('toimii')
    else:
        mapping_dict = {}
        filepath = filename.path
        linecount = 0
        with open(filepath, 'r') as f:
            for line in f:
                linecount +=1  
                try:
                    generic_name = line.split()[0]
                except IndexError:
                    print(filename.path)
                    exit()
                        
                mapping_dict[generic_name] = {}
                try:
                    mapping_dict[generic_name]['ATOMNAME'] = line.split()[1]
                except IndexError:
                    print(filename.path)
                    exit()
                try:
                    mapping_dict[generic_name]['RESIDUE'] = line.split()[2]
                except IndexError:
                    pass
                    
                # PHOSPHOLIPIDS: automatically assign atoms to headgroup, glycerol backbone, sn-1 or sn-2 
                # CHOLESTEROL: headgroup is OH and the rest is tail
                # does not work for other molecules 
                fragment = ""
                    
                if 'CHOL' in filename.path:
                    if 'M_C3O' in generic_name:
                        fragment = 'headgroup'
                    else:
                        fragment = 'tail'
                elif 'PC' in filename.path or 'PE' in filename.path or 'PS' in filename.path or 'PG' in filename.path or 'PI' in filename.path: #phospholipids
                    if 'M_G3_M' in generic_name or 'M_G3H' in generic_name or 'M_G1_M' in generic_name or 'M_G1H' in generic_name or 'M_G2_M' in generic_name or 'M_G2H' in generic_name: 
                        fragment = 'glycerol backbone'
                    elif 'M_G3' in generic_name:
                        fragment = 'headgroup'
                    elif 'M_G1C' in generic_name:
                        fragment = 'sn-1'
                    elif 'M_G2C' in generic_name:
                        fragment = 'sn-2'
                else:
                    fragment = ""
                        
                if fragment:
                    mapping_dict[generic_name]['FRAGMENT'] = fragment

            f.close()
        
        #save mapping file as a dictionary in yaml format
        fname = filename.path.replace('.txt','').split('/')[-1]
        # check that all atoms are included
        print("number of atoms : " + str(len(mapping_dict.keys())) + "         " + fname)
        print("number of atoms in txt file: " + str(linecount))
        print("####################")
        
        if len(mapping_dict.keys()) == linecount:
            new_file = './mapping_files/mapping_files_yaml/' + fname + '.yaml'
           # print(new_file)
        
            with open(new_file,'w') as f:
                yaml.dump(mapping_dict, f)
            f.close()

