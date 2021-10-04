import os
import sys
import numpy as np
import json
import matplotlib.pyplot as plt
import MDAnalysis
import urllib.request
import yaml

sys.path.insert(1, '../BuildDatabank/')
from databankLibrary import download_link, lipids_dict, molecule_numbers_dict,databank, read_trajs_calc_OPs, parse_op_input, find_OP, OrderParameter
import buildH_calcOP_test

path = '../../Data/Simulations/'
db_data = databank(path)
systems = db_data.get_systems()

ready = 0
notready = 0
for system in systems:
    Nlipid = 0
    path = system['path']

    for key in system['COMPOSITION']:
        outfilename = path + key + 'OrderParameters.json'
        #print(outfilename)
        if os.path.isfile(outfilename):
            FileFound = True
        elif key in lipids_dict:
            FileFound = False
            continue

    if FileFound:
        continue

    print('Analyzing: ', system['path'])

    # This is temprorary because not all systems are converted to new format yet
    # To be removed
    try:    
        for molecule in system['COMPOSITION']:
            if molecule in lipids_dict:
                Nlipid += np.sum(system['COMPOSITION'][molecule]['COUNT'])
    except:
        notready = notready + 1
        print(system['path'], ' not ready')
        continue

                
    doi = system.get('DOI')
    trj = system.get('TRJ')
    tpr = system.get('TPR')
    trj_name = path + system.get('TRJ')[0][0]
    tpr_name = path + system.get('TPR')[0][0]
    trj_url = download_link(doi, trj[0][0])
    tpr_url = download_link(doi, tpr[0][0])
    #print(trj_name,tpr_name)
    
    if (not os.path.isfile(tpr_name)):
        response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
    if (not os.path.isfile(trj_name)):
        response = urllib.request.urlretrieve(trj_url, trj_name)

    software=system['SOFTWARE']
    EQtime=float(system['TIMELEFTOUT'])*1000
    try:
        unitedAtom = system['UNITEDATOM_DICT']
    except:
        unitedAtom = False

    if 'gromacs' in software:
         xtcwhole= path + '/whole.xtc'
         print("Make molecules whole in the trajectory")
         os.system('echo System | gmx trjconv -f ' + trj_name + ' -s ' + tpr_name + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime))
    else:
        print('Order parameter calculation for other gromacs is yet to be implemented.')
        continue


    if unitedAtom:
        for key in system['UNITEDATOM_DICT']:
        #construct order parameter definition file for CH bonds from mapping file
            def_file = open(path + key + '.def', 'w')

            mapping_file = system['COMPOSITION'][key]['MAPPING']
            previous_line = ""
            
            with open('../BuildDatabank/mapping_files/'+mapping_file, "r") as f:
                for line in f.readlines():
                    if not line.startswith("#"):
                        regexp1_H = re.compile(r'M_[A-Z0-9]*C[0-9]*H[0-9]*_M')
                        regexp2_H = re.compile(r'M_G[0-9]*H[0-9]*_M')
                        regexp1_C = re.compile(r'M_[A-Z0-9]*C[0-9]*_M')
                        regexp2_C = re.compile(r'M_G[0-9]_M')

                        if regexp1_C.search(line) or regexp2_C.search(line):
                            atomC = line.split()
                            atomH = []
                        elif regexp1_H.search(line) or regexp2_H.search(line):
                            atomH = line.split()
                        else:
                            atomC = []
                            atomH = []

                        if atomH:
                            items = [atomC[1], atomH[1], atomC[0], atomH[0]]
                            def_line = items[2] + " " + items[3] + " " + key + " " + items[0] + " " + items[1] + "\n"
                            if def_line != previous_line:
                                def_file.write(def_line)
                                print(def_line)
                                previous_line = def_line
            def_file.close()
        #Add hydrogens to trajectory and calculate order parameters with buildH
        ordPfile = path + key + 'OrderParameters.dat' 
        topfile = path + '/frame0.gro'
                        
        os.system('echo System | gmx trjconv -f ' + xtcwhole + ' -s ' + tpr_name + ' -dump 0 -o ' + topfile )
                         
        deffile = path + key + '.def' 
        lipidname = system['UNITEDATOM_DICT'][key]
        #    print(lipidname)
        buildH_calcOP_test.main(topfile,lipidname,deffile,xtcwhole,ordPfile)

        outfile=open(ordPfile,'w')
        line1="Atom     Average OP     OP stem"+'\n'
        outfile.write(line1)
        
        data = {}
        outfile2= path + key + 'OrderParameters.json'
        
        with open(ordPfile + '.jmelcr_style.out') as OPfile:
            lines=OPfile.readlines()
            for line in lines:
                if "#" in line:
                    continue
                line2 = line.split()[0] + " " + line.split()[1] + "  " + line.split()[5] + "  " + line.split()[6] + " " + line.split()[7] + "\n"
                outfile.write(line2)

                OPname = line.split()[0] + " " + line.split()[1]
                OPvalues = [line.split()[5], line.split()[6] ,line.split()[7]]
                data[str(OPname)]=[]
                data[str(OPname)].append(OPvalues)
        
            with open(outfile2, 'w') as f:
                json.dump(data,f)

        outfile.close()
        outfile.close()

        # os.system('cp ' + str(dir_tmp) + '/' + key + 'OrderParameters.dat ' + DATAdir) #Or should these be put into Data/Simulations/
        # os.system('cp ' +str(dir_tmp) + '/' + key + 'OrderParameters.json ' + DATAdir)
    else:
        #trj = str(DATAdir) + '/' + str(trj)
        gro = path + '/conf.gro'

        #make gro file
        print("\n Makin gro file")
        os.system('echo System | gmx trjconv -f ' + trj_name + ' -s ' + tpr_name + ' -dump 0 -o ' + gro)
                    
        for key in system['COMPOSITION']:
            if key in lipids_dict.keys(): 
                mapping_file = system['COMPOSITION'][key]['MAPPING']
                resname = system['COMPOSITION'][key]['NAME']
                OrdParam=find_OP('../BuildDatabank/mapping_files/'+mapping_file,gro,xtcwhole,resname)

                outfile=open(path + key + 'OrderParameters.dat','w')
                line1="Atom     Average OP     OP stem"+'\n'
                outfile.write(line1)
    
                data = {}
                outfile2=path + key + 'OrderParameters.json' 

                for i,op in enumerate(OrdParam):
                    resops =op.get_op_res
                    (op.avg, op.std, op.stem) =op.get_avg_std_stem_OP
                    line2=str(op.name)+" "+str(op.avg)+" "+str(op.stem)+'\n'
                    outfile.write(line2)
    
                    data[str(op.name)]=[]
                    data[str(op.name)].append(op.get_avg_std_stem_OP)
        
                with open(outfile2, 'w') as f:
                    json.dump(data,f)
                outfile.close()
                f.close()
                # os.system('cp ' + str(dir_path) + '/' + key + 'OrderParameters.dat ' + DATAdir) #MUUTA
                #os.system('cp ' +str(dir_path) + '/' + key + 'OrderParameters.json ' + DATAdir) #MUUTA
    
    print("Order parameters calculated and saved to " + DATAdir)

    ready = ready + 1
        
print('ready',ready)
print('notready',notready)
