import os
import sys
import yaml
import json
import matplotlib.pyplot as plt
import numpy as np
import math
import argparse
import decimal as dc

from random import randint

from matplotlib import cm
from scipy.stats import norm

import urllib.request
from urllib.error import URLError,HTTPError,ContentTooShortError

sys.path.insert(1, '../BuildDatabank/')
from databankLibrary import download_link, lipids_dict, read_trajs_calc_OPs, parse_op_input, find_OP, OrderParameter
import buildH_calcOP_test


#order parameters or quality analysis
parser = argparse.ArgumentParser(description="")
parser.add_argument("-op",action='store_true', help="Calculate order parameters for simulations in data bank "
                        "format.")
parser.add_argument("-q",action='store_true', help="Quality evaluation of simulated order parameters "
                        "format.")
args = parser.parse_args()

lipid_numbers_list = lipids_dict.keys() # should contain all lipid names
#################################
class Simulation:
    def __init__(self, readme, data, indexingPath):
        self.readme = readme
        self.data = data #dictionary where key is the lipid type and value is order parameter file
        self.indexingPath = indexingPath
        
    def getLipids(self, molecules=lipid_numbers_list):
        lipids = []
        for key in molecules:
            try:
                if self.readme['N'+key] != [0,0]: 
                    lipids.append(key)
            except KeyError:
                continue
        return lipids     
        
class Experiment:
    pass
    

#Quality evaluation of simulated data
#assume that an order parameter calculated S from a simulation is normally distributed
# OP_sd = sqrt(N)*STEM, where N is number of lipids and STEM is standard error of mean
#def OP_sd(N, STEM):
#    op_sd = math.sqrt(N)*STEM
#    
#    return op_sd

# op_sd = STEM
    
# P: what is the probability that S_exp +/- 0.02 is in g(s) where g(s) is the probability density function of normal distribution N(s, S_sim, S_sim_sd)

def prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_sd):
    #normal distribution N(s, OP_sim, op_sim_sd)
    a = OP_exp - exp_error
    b = OP_exp + exp_error
    #P_S = norm.cdf(b, loc=OP_sim, scale=op_sim_sd) - norm.cdf(a, loc=OP_sim, scale=op_sim_sd)
    #changing to survival function to increase precision, note scaling. Must be recoded to increase precision still
    P_S = -norm.sf(b, loc=OP_sim, scale=op_sim_sd) + norm.sf(a, loc=OP_sim, scale=op_sim_sd)

    if math.isnan(P_S) :
     return P_S

    #this is an attempt to deal with precision, max set manually to 70
    dc.getcontext().prec=70
    precise_log=-dc.Decimal(P_S).log10()

    #print("difference of two prec logs",float(foo)+math.log10(P_S))

    return float(precise_log)

# quality of simulated order parameter
#def OPquality(P_S,op_sim_STEM):
#   # print('probability')
#    #print(P_S)
#    if P_S != 0:
#        quality = 1-math.log(P_S) #/op_sim_STEM     # math.log(P_S/op_sim_STEM)      #/ math.sqrt(op_sim_STEM) #/ (op_sim_STEM*op_sim_STEM)
#    else:
#        quality = 0
   # quality_float = quality.item()
 #   print('quality')
 #   print(quality)
#    return quality
    
#
def OPquality(OP_exp, OP_sim):
    quality = np.absolute(OP_exp - OP_sim)
    return quality
    
# quality of molecule fragments

def evaluated_percentage(fragment, exp_op_data):
    count_value = 0
    fragment_size = 0
    for key, value in exp_op_data.items():
        for f in fragment:
            if f in key:
                fragment_size += 1
                if value[0][0] != 'nan':
                    count_value += 1
    if fragment_size != 0:
        return count_value / fragment_size
    else:
        return 0
    
def carbonError(OP_sim, OP_exp):
    
    E_i = 0
    quality = OPquality(OP_exp, OP_sim)

    if quality > 0.02:
        E_i = quality - 0.02
    else:
        E_i = 0
    return E_i  



def fragmentQuality(fragment, exp_op_data, sim_op_data):
    p_F = evaluated_percentage(fragment, exp_op_data)
    #warning, hard coded experimental error
    exp_error=0.02
    E_sum = 0
    if p_F != 0:
        for key_exp, value_exp in exp_op_data.items():
            #  print(key_exp)
            #  print(value_exp)
            if fragment in key_exp and value_exp[0][0] != 'nan':
                OP_exp = value_exp[0][0]
               # print(OP_exp)
                OP_sim = sim_op_data[key_exp][0]
               # print(OP_sim)
               # op_sim_STEM=sim_op_data[key_exp][0][2]
                E_sum += carbonError(OP_exp, OP_sim)
               #change here if you want to use shitness(TM) scale for fragments. Warning big umbers will dominate
               # E_sum+= prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_STEM)
        E_F = E_sum / p_F
        return E_F
    else:
        return 'nan'

    
###################################################################################################
simulations = []
for subdir, dirs, files in os.walk(r'../../Data/Simulations/'): #
    for filename1 in files:
        filepath = subdir + os.sep + filename1
        
        if filepath.endswith("README.yaml"):
            READMEfilepathSimulation = subdir + '/README.yaml'
            with open(READMEfilepathSimulation) as yaml_file_sim:
                readmeSim = yaml.load(yaml_file_sim, Loader=yaml.FullLoader)
                indexingPath = "/".join(filepath.split("/")[4:8])
                print(indexingPath)
                print(filepath)
                try:
                    if readmeSim['EXPERIMENT']:

                        simOPdata = {} #order parameter files for each type of lipid
                        for filename2 in files:
                            if filename2.endswith('OrderParameters.json'):
                                key_data1 = filename2.replace('OrderParameters.json', '')
                                OPfilepath = subdir + "/" + filename2
                                with open(OPfilepath) as json_file:
                                    simOPdata[key_data1] = json.load(json_file)
                                    json_file.close()
                        simulations.append(Simulation(readmeSim, simOPdata, indexingPath))
                        yaml_file_sim.close()
                except KeyError:
                    # print("No matching experimental data for system " + readmeSim['SYSTEM'] + " in directory " + indexingPath)
                    continue
                    

if (not os.path.isdir('../../Data/QualityEvaluation/')): 
    os.system('mkdir ../../Data/QualityEvaluation/')

for simulation in simulations:
    sub_dirs = simulation.indexingPath.split("/")
    os.system('mkdir ../../Data/QualityEvaluation/' + sub_dirs[0])
    os.system('mkdir ../../Data/QualityEvaluation/' + sub_dirs[0] + '/' + sub_dirs[1])
    os.system('mkdir ../../Data/QualityEvaluation/' + sub_dirs[0] + '/' + sub_dirs[1] + '/' + sub_dirs[2])    
    os.system('mkdir ../../Data/QualityEvaluation/' + sub_dirs[0] + '/' + sub_dirs[1] + '/' + sub_dirs[2] + '/' + sub_dirs[3])
    
    #save fragment quality here
    DATAdir_FQ = '../../Data/QualityEvaluation/' + str(sub_dirs[0]) + '/' + str(sub_dirs[1]) + '/' + str(sub_dirs[2]) + '/' + str(sub_dirs[3])
    
    #save OP quality here
    DATAdir_OP = '../../Data/Simulations/' + + str(sub_dirs[0]) + '/' + str(sub_dirs[1]) + '/' + str(sub_dirs[2]) + '/' + str(sub_dirs[3])
   # print(DATAdir)
     

    for lipid1 in simulation.getLipids():
        #print(lipid1)
        print(simulation.indexingPath)
        #print(simulation.data.keys())
        
        # OP_data_lipid = simulation.data[lipid1]
        OP_data_lipid = {}
        #convert elements to float because in some files the elements are strings
        for key, value in simulation.data[lipid1].items():
            OP_array = [float(x) for x in simulation.data[lipid1][key][0]]  
            OP_data_lipid[key] = OP_array
            
        OP_qual_data = {}
        # go through file paths in simulation.readme['EXPERIMENT']
        print(simulation.readme['EXPERIMENT'].values())
        for value in simulation.readme['EXPERIMENT'].values():
          # get readme file of the experiment
            experimentFilepath = "../../Data/experiments/" + value
            print(experimentFilepath)
            READMEfilepathExperiment  = experimentFilepath + '/README.yaml'
            experiment = Experiment()
            with open(READMEfilepathExperiment) as yaml_file_exp:
                readmeExp = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)
                experiment.readme = readmeExp
                #print(experiment.readme)
            yaml_file_exp.close()

            lipidExpOPdata = {}
            try:
                exp_OP_filepath = experimentFilepath + '/' + lipid1 + '_Order_Parameters.json'
            except FileNotFoundError:
                print("Experimental order parameter data do not exist for lipid " + lipid1 + ".")
                continue
            else:
                #print(exp_OP_filepath)
                with open(exp_OP_filepath) as json_file:
                    lipidExpOPdata = json.load(json_file)
                json_file.close()
            
                simulationREADMEsave = DATAdir + '/README.yaml'
                with open(simulationREADMEsave, 'w') as f:
                    yaml.dump(simulation.readme,f, sort_keys=False)
                f.close()

            exp_error = 0.02
            
            for key in OP_data_lipid.keys():
                OP_array = OP_data_lipid[key]
                try:
                    OP_exp = lipidExpOPdata[key][0][0]
                except KeyError:
                    OP_array.append('NaN')
                    continue
                else:
                    if lipidExpOPdata[key][0][0] is not 'NaN':
                        OP_sim = OP_array[0]
                        op_sim_STEM = OP_array[2]
                        #changing to use shitness(TM) scale. This code needs to be cleaned
                        op_quality = prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_STEM)
                        OP_array.append(op_quality)
                    else:
                        OP_array.append('NaN')
                
                OP_qual_data[key] = OP_array    
                
             print(OP_qual_data)
#            for key, value in lipidExpOPdata.items():
#                if lipidExpOPdata[key][0][0] is not 'NaN':
#                    OP_array = OP_data_lipid[key] #[float(x) for x in OP_data_lipid[key][0]] #convert elements to float because in some files the elements are strings 
#                    print(OP_array)
#                    #print(type(OP_array))
#                    OP_exp = value[0][0]
#                    OP_sim = OP_array[0]
#                    #op_sim_sd = OP_array[1] 
#                    op_sim_STEM = OP_array[2] 
#                    #changing to use shitness(TM) scale. This code needs to be cleaned
#                    op_quality = prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_STEM) #(OP_exp, exp_error, OP_sim, op_sim_sd)
#             
#                    # op_quality = OPquality(OP_exp, OP_sim)   #numpy float must be converted to float
#                    # print(type(op_quality))
#                    OP_array.append(op_quality)
#                    #print(OP_array)
#                
#                    OP_qual_data[key] = OP_array

#            print(OP_qual_data) 
        # quality data should be written into the OrderParameters.json file of the simulation                  
            outfile = DATAdir_OP + '/' + lipid1 + '_OrderParameters.json'
        
            with open(outfile, 'w') as f:
                json.dump(OP_qual_data,f)
            f.close()
            
        # calculate quality for molecule fragments headgroup, sn-1, sn-2
            headgroup = fragmentQuality(['M_G3','M_G1_','M_G2_'], lipidExpOPdata, OP_data_lipid)
            sn1 = fragmentQuality(['M_G1C'], lipidExpOPdata, OP_data_lipid)
            sn2 = fragmentQuality(['M_G2C'], lipidExpOPdata, OP_data_lipid)
            
            fragment_quality = {}
            fragment_quality['headgroup'] = headgroup
            fragment_quality['sn-1'] = sn1
            fragment_quality['sn-2'] = sn2
           
            print('headgroup ')
            print(headgroup)
            print('sn1 ') 
            print(sn1)
            print('sn2 ') 
            print(sn2) 
            
            fragment_quality_file = DATAdir_FQ + '/' + lipid1 + '_FragmentQuality.json'
            
            with open(fragment_quality_file, 'w') as f:
                json.dump(fragment_quality,f)
            f.close()
        
        
        
     #   print(OP_qual_data)                        
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
