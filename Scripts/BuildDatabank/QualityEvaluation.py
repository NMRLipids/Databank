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

lipid_numbers_list = lipids_dict.keys() # should contain all lipid names
#################################
class Simulation:
    def __init__(self, readme, OPdata, FFdata,indexingPath):
        self.readme = readme
        self.OPdata = OPdata #dictionary where key is the lipid type and value is order parameter file
        self.FFdata = FFdata
        self.indexingPath = indexingPath
        
    def getLipids(self, molecules=lipid_numbers_list):
        lipids = []
        
        for key in self.readme['COMPOSITION'].keys():
            if key in molecules:
                lipids.append(key)
        return lipids
        
    def molarFraction(self, molecule,molecules=lipid_numbers_list): #only for lipids
        sum_lipids = 0
        number = sum(self.readme['COMPOSITION'][molecule]['COUNT']) 
        
        for key in self.readme['COMPOSITION'].keys():
            if key in molecules:
                sum_lipids += sum(self.readme['COMPOSITION'][key]['COUNT'])

        return number / sum_lipids
        
class Experiment:
    pass


#Quality evaluation of simulated data
#Order parameters
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
    
# quality of molecule fragments

def evaluated_percentage(fragments, exp_op_data):
    count_value = 0
    fragment_size = 0
    for key, value in exp_op_data.items():
        for f in fragments:
             if f in key:
                 fragment_size += 1
                 if value[0][0] != 'nan':
                     count_value += 1
    if fragment_size != 0:
        return count_value / fragment_size
    else:
        return 0


def fragmentQuality(fragments, exp_op_data, sim_op_data):
    #print("fragment quality")
    p_F = evaluated_percentage(fragments, exp_op_data)
    #warning, hard coded experimental error
    exp_error=0.02
    E_sum = 0
    AV_sum = 0
    if p_F != 0:
        for fr in fragments:
            for key_exp, value_exp in exp_op_data.items():
            #    print(key_exp)
            #    print(value_exp[0][0])
                if (fr in key_exp) and value_exp[0][0] != 'nan':
                    OP_exp = value_exp[0][0]
                    # print(OP_exp)
                    OP_sim = sim_op_data[key_exp][0]
                    # print(OP_sim)

                    #print(sim_op_data[key_exp])
                    op_sim_STEM=sim_op_data[key_exp][2]
                    
                    #change here if you want to use shitness(TM) scale for fragments. Warning big numbers will dominate
                    #if OP_exp != 'NaN':
                    QE = prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_STEM)
                    #print(prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_STEM))
                    if QE >0: # 0! = 'nan': #QE > 0 and  QE != 'inf': # and  QE != 'nan'
                        #if QE == 'NaN':
                        #    E_sum = E_sum
                        if QE == float("inf"): #'Infinity' or QE == 'inf':
                            E_sum += 300
                            AV_sum += 1
                        else:
                            #print(QE)
                            E_sum += prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_STEM)
                            AV_sum += 1

        E_F = (E_sum / AV_sum) / p_F
        return E_F
    else:
        return 'nan'
        
def fragmentQualityAvg(lipid,fragment_qual_dict):
    if lipid != 'CHOL':
        headgroup_sum = 0
        sn1_sum = 0
        sn2_sum = 0

        headgroup_c = 0
        sn1_c = 0
        sn2_c = 0

        for doi in fragment_qual_dict.keys():
            #print(doi)
            for key in fragment_qual_dict[doi].keys():
                if key == 'headgroup' and fragment_qual_dict[doi][key] != 'nan':
                    headgroup_sum += fragment_qual_dict[doi][key]
                    headgroup_c += 1
                elif key == 'sn-1' and fragment_qual_dict[doi][key] != 'nan':
                    sn1_sum += fragment_qual_dict[doi][key]
                    sn1_c += 1
                elif key == 'sn-2' and fragment_qual_dict[doi][key] != 'nan':
                    sn2_sum += fragment_qual_dict[doi][key]
                    sn2_c += 1

        if headgroup_sum != 0:
            headgroup_avg = headgroup_sum / headgroup_c
        else:
            headgroup_avg = 'nan'

        if sn1_sum != 0:
            sn1_avg = sn1_sum / sn1_c
        else:
            sn1_avg = 'nan'

        if sn2_sum != 0:
            sn2_avg = sn2_sum / sn2_c
        else:
            sn2_avg = 'nan'

        total_quality = 0
        if headgroup_avg != 'nan' and sn1_avg != 'nan' and sn2_avg != 'nan':
            total_quality = (headgroup_avg + sn1_avg + sn2_avg) / 3
        else:
            total_quality = 'nan'

        return  headgroup_avg, sn1_avg, sn2_avg, total_quality

    else:
        qual_sum = 0
        chol_c = 0
        for doi in fragment_qual_dict.keys():
            qual_sum += fragment_qual_dict[doi]['cholesterol']
            chol_c += 1
        if chol_c > 0:
            total_quality = qual_sum / chol_c
        else:
            total_quality = 0
        return total_quality
        

def systemQuality(system_fragment_qualities):
    out_dict = {}
    headgroup = []
    sn1 = []
    sn2 = []
    total = []
     
    w_nan = []

    for lipid in system_fragment_qualities.keys():
        w = simulation.molarFraction(lipid)
        if lipid != 'CHOL':    # SHOULD BE CHANGED TO WORK ALSO WITH OTHER LIPIDS WITHOUT HEAD AND TAILS THAN CHOLESTEROL
            for key, value in system_fragment_qualities[lipid].items():
                if value != 'nan':
                    if key == 'headgroup':
                        headgroup.append(w * value)
                    elif key == 'sn-1':
                        sn1.append(w * value)
                    elif key == 'sn-2':
                        sn2.append(w * value)
                    elif key == 'total':
                        total.append(w * value)
                    else:
                        continue
                else:
                    w_nan.append(1-w) # save 1 - w of a lipid into a list if the fragment quality is nan
        else:
            for key, value in system_fragment_qualities[lipid].items():
                 if value != 'nan':
                     if key == 'total':
                         total.append(w * value)

    out_dict['headgroup'] = sum(headgroup) / np.prod(w_nan) # multiply all elements of w_nan and divide the sum by the product
    out_dict['sn-1'] = sum(sn1) / np.prod(w_nan)
    out_dict['sn-2'] = sum(sn2) / np.prod(w_nan)
    out_dict['total'] = sum(total) / np.prod(w_nan)

    ## EXTREMELY DIRTY FIX FOR WORKSHOP, SHOULD BE IMPROVED LATER
#    for lipid in system_fragment_qualities.keys():
#        w = simulation.molarFraction(lipid)
#        if lipid != 'CHOL':    
#            for key, value in system_fragment_qualities[lipid].items():                        
#                if value == 'nan':
#                   if key == 'headgroup':
#                       headgroup[:] = [x / (1-w) for x in headgroup]
#                   elif key == 'sn-1':
#                       sn1[:] = [x / (1-w) for x in sn1] 
#                   elif key == 'sn-2':
#                       sn2[:] = [x / (1-w) for x in sn2] 
#                   elif key == 'total':
#                       total[:] = [x / (1-w) for x in total]
#                else:
#                    continue

                         
#    out_dict['headgroup'] = sum(headgroup)
#    out_dict['sn-1'] = sum(sn1)
#    out_dict['sn-2'] = sum(sn2)
#    out_dict['total'] = sum(total)

    return out_dict
    
#Form factor quality


def calc_k_e(simFFdata,expFFdata):
    """Scaling factor as defined by Kučerka et al. 2008b, doi:10.1529/biophysj.107.122465  """
    sum1 = 0
    sum2 = 0
    
   # print("simulation:" + str(len(simFFdata)))
   # print("experiment:" + str(len(expFFdata)))
   
    if len(expFFdata) <= len(simFFdata):
        for i in range(0,len(expFFdata)): #experiment should contain less data points
            F_s = simFFdata[i][1]
            F_e = expFFdata[i][1]
            deltaF_e = expFFdata[i][2]
        
            sum1 = sum1 + np.abs(F_s)*np.abs(F_e)/(deltaF_e**2)
            sum2 = sum2 + np.abs(F_e)**2 / deltaF_e**2
        k_e = sum1 / sum2
        return k_e
    
    else:
        return ""
    
    


def formfactorQuality(simFFdata, expFFdata):
    """Calculate form factor quality for a simulation as defined by Kučerka et al. 2010, doi:10.1007/s00232-010-9254-5 """
    k_e = calc_k_e(simFFdata,expFFdata)
    N = len(expFFdata)
    
    sum1 = 0
    
    if len(expFFdata) <= len(simFFdata):
        for i in range(0,len(expFFdata)): #experiment should contain less data points
            F_s = simFFdata[i][1]
            F_e = expFFdata[i][1]
            deltaF_e = expFFdata[i][2] 
        
            sum1 = sum1 + (np.abs(F_s) - k_e*np.abs(F_e))**2 / deltaF_e**2
    
        khi2 = np.sqrt(sum1) / np.sqrt(N - 1)
    
        return khi2
    else:
        
        return ""

      

def loadSimulations():
    simulations = []
    for subdir, dirs, files in os.walk(r'../../Data/Simulations/'): #
        for filename1 in files:
            filepath = subdir + os.sep + filename1
        
            if filepath.endswith("README.yaml"):
                READMEfilepathSimulation = subdir + '/README.yaml'
                readmeSim = {}
                with open(READMEfilepathSimulation) as yaml_file_sim:
                    readmeSim = yaml.load(yaml_file_sim, Loader=yaml.FullLoader)
                yaml_file_sim.close()    
                indexingPath = "/".join(filepath.split("/")[4:8])

                #print(indexingPath)
                #print(filepath)
                #print(readmeSim)

                try:
                    experiments = readmeSim['EXPERIMENT']
                except KeyError:
                    # print("No matching experimental data for system " + readmeSim['SYSTEM'] + " in directory " + indexingPath)
                    continue
                else:
                    if any(experiments.values()): #if experiments is not empty
                        #print(any(experiments))
                        #print('Experiments found for' + filepath)
                        #print(experiments)
                        simOPdata = {} #order parameter files for each type of lipid
                        simFFdata = {} # form factor data
                        for filename2 in files:
                            if filename2.endswith('OrderParameters.json'):
                                lipid_name = filename2.replace('OrderParameters.json', '')
                             #  print(lipid_name)
                                dataPath = subdir + "/" + filename2
                                OPdata = {}
                                with open(dataPath) as json_file:
                                    OPdata = json.load(json_file)
                                json_file.close()
                                simOPdata[lipid_name] = OPdata
                                
                            elif filename2 == "FormFactor.json":
                                dataPath = subdir + "/" + filename2
                                with open(dataPath) as json_file:
                                    simFFdata = json.load(json_file)
                                json_file.close()
                                    
                        simulations.append(Simulation(readmeSim, simOPdata, simFFdata, indexingPath))
                    else:
                        #print("The simulation does not have experimental data.")
                        continue
                
                    
    return simulations



###################################################################################################
print('start')
simulations = loadSimulations()

#if (not os.path.isdir('../../Data/QualityEvaluation/')): 
#    os.system('mkdir ../../Data/QualityEvaluation/')


for simulation in simulations:
    sub_dirs = simulation.indexingPath.split("/")
    
    #save OP quality and FF quality here
    DATAdir = '../../Data/Simulations/' + str(sub_dirs[0]) + '/' + str(sub_dirs[1]) + '/' + str(sub_dirs[2]) + '/' + str(sub_dirs[3])
    #print(DATAdir)
   
   
    #Order Parameters 
    system_quality = {}
    for lipid1 in simulation.getLipids():
        #print(lipid1)
        print('Simulation path ' + simulation.indexingPath)
        #print(simulation.OPdata.keys())
        
        # OP_data_lipid = simulation.OPdata[lipid1]
        OP_data_lipid = {}
        #convert elements to float because in some files the elements are strings
        for key, value in simulation.OPdata[lipid1].items():
            OP_array = [float(x) for x in simulation.OPdata[lipid1][key][0]]  
            OP_data_lipid[key] = OP_array
            
        
        
        # go through file paths in simulation.readme['EXPERIMENT']
        #print(simulation.readme['EXPERIMENT'].values())

        
        for lipid, experiments in simulation.readme['EXPERIMENT']['ORDERPARAMETER'].items():
            print(lipid,experiments)
            data_dict = {}
            fragment_qual_dict = {}
            for doi, path in experiments.items():
                OP_qual_data = {}
            # get readme file of the experiment
                experimentFilepath = "../../Data/experiments/OrderParameters/" + path
                print('Experimental path ' + experimentFilepath)
                READMEfilepathExperiment  = experimentFilepath + '/README.yaml'
                experiment = Experiment()
                with open(READMEfilepathExperiment) as yaml_file_exp:
                    readmeExp = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)
                    experiment.readme = readmeExp
                    #print(experiment.readme)
                yaml_file_exp.close()

                exp_OP_filepath = experimentFilepath + '/' + lipid1 + '_Order_Parameters.json'
                #print(exp_OP_filepath)
                lipidExpOPdata = {}
                try:
                    with open(exp_OP_filepath) as json_file:
                        lipidExpOPdata = json.load(json_file)
                    json_file.close()
                except FileNotFoundError:
                    print("Experimental order parameter data do not exist for lipid " + lipid1 + ".")
                    break


                exp_error = 0.02
           
                for key in OP_data_lipid.keys():
                    OP_array = OP_data_lipid[key].copy()
                    try:
                        OP_exp = lipidExpOPdata[key][0][0]
                    except KeyError:
     #                   OP_array.append('NaN')
     #                   OP_qual_data[key] = OP_array
                        continue
                    else:
                        if OP_exp is not 'NaN':
                            OP_sim = OP_array[0]
                            op_sim_STEM = OP_array[2]
                            #changing to use shitness(TM) scale. This code needs to be cleaned
                            op_quality = prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_STEM)
                            OP_array.append(OP_exp)
                            OP_array.append(exp_error)   #hardcoded!!!! 0.02 for all experiments
                            OP_array.append(op_quality)
                      #  else:
                        #    OP_array.append(OP_exp)
                        #    OP_array.append(exp_error)  #hardcoded!!!! 0.02 for all experiments
                       #     OP_array.append('NaN')
                
                    OP_qual_data[key] = OP_array    
                
                #print(OP_qual_data)
                
                data_dict[doi] = OP_qual_data
                
                # calculate quality for molecule fragments headgroup, sn-1, sn-2

                fragment_quality = {}

                if lipid1 == 'CHOL':
                    cholQ = fragmentQuality(['M_C'], lipidExpOPdata, OP_data_lipid)
                    fragment_quality['cholesterol'] = cholQ
                else:
                    headgroup = fragmentQuality(['M_G3','M_G1_','M_G2_'], lipidExpOPdata, OP_data_lipid)
                    sn1 = fragmentQuality(['M_G1C'], lipidExpOPdata, OP_data_lipid)
                    sn2 = fragmentQuality(['M_G2C'], lipidExpOPdata, OP_data_lipid)

                    fragment_quality['headgroup'] = headgroup
                    fragment_quality['sn-1'] = sn1
                    fragment_quality['sn-2'] = sn2

                fragment_qual_dict[doi] = fragment_quality
                
            fragment_quality_output = {}
            if lipid1 != 'CHOL':
                headgroup_avg, sn1_avg, sn2_avg, total_qual = fragmentQualityAvg(lipid1,fragment_qual_dict)
                fragment_quality_output['headgroup'] = headgroup_avg
                fragment_quality_output['sn-1'] = sn1_avg
                fragment_quality_output['sn-2'] = sn2_avg
                fragment_quality_output['total'] = total_qual
            else:
               # print("Cholesterol works")
                total_qual = fragmentQualityAvg(lipid1,fragment_qual_dict)
                fragment_quality_output['total'] = total_qual
            
            system_quality[lipid1] = fragment_quality_output

            fragment_quality_file = DATAdir + '/' + lipid1 + '_FragmentQuality.json'

            FGout = False
            for FG in fragment_quality_output:
                #print(FG,fragment_quality_output[FG])
                if fragment_quality_output[FG] == 'nan':
                    continue
                if fragment_quality_output[FG] > 0:
                    FGout = True
            if FGout:
                with open(fragment_quality_file, 'w') as f:
                    json.dump(fragment_quality_output,f)
                f.close()

                
                
                
        

            #write into the OrderParameters_quality.json quality data file                  
            outfile1 = DATAdir + '/' + lipid1 + '_OrderParameters_quality.json'
            #doi : {'carbon hydrogen': [op_sim, sd_sim, stem_sim, op_exp, exp_error, quality] ... }
            if(len(data_dict) > 0):
                with open(outfile1, 'w') as f:
                    json.dump(data_dict,f)
                f.close()

        print('input to system quality')
        print(system_quality)
        #calculate system quality
        system_qual_output = systemQuality(system_quality)
        print('system')
        print(system_qual_output)
        #make system quality file
        outfile2 = DATAdir + '/SYSTEM_quality_test.json'
        SQout = False
        for SQ in system_qual_output:
            if system_qual_output[SQ] > 0:
                SQout = True
        if SQout:
            with open(outfile2, 'w') as f:
                json.dump(system_qual_output,f)
            f.close() 
        
        print('')
        
        
     #   print(OP_qual_data)                        
                
                
    #Form factor quality
        
    expFFpath = simulation.readme['EXPERIMENT']['FORMFACTOR']
    #print(ffexpPath)
    expFFdata = {}
    for subdir, dirs, files in os.walk(r'../../Data/experiments/FormFactors/' + expFFpath + '/'):
        for filename in files:
            filepath = '../../Data/experiments/FormFactors/' + expFFpath + '/' + filename
            if filename.endswith('_FormFactor.json'):
                print(filename)
                with open(filepath) as json_file:
                    expFFdata = json.load(json_file)
                json_file.close()
    
    
    simFFdata = simulation.FFdata
    
    ffQuality = formfactorQuality(simFFdata, expFFdata)
    
    outfile3 = DATAdir + '/FormFactorQuality.json'
    with open(outfile3,'w') as f:
        json.dump(ffQuality,f)
    f.close()
    
          
                
                
                
                
                
                
                
                
                
                
                
                
                
                
