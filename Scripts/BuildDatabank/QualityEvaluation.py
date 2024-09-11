#!/usr/bin/env python3
# coding: utf-8

import os, sys
import yaml, json
import numpy as np


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from DatabankLib import NMLDB_SIMU_PATH, NMLDB_EXP_PATH
from DatabankLib.jsonEncoders import CompactJSONEncoder
import DatabankLib.quality as qq

#################################

simulations = qq.loadSimulations()

EvaluatedOPs = 0
EvaluatedFFs = 0

for simulation in simulations:
    #save OP quality and FF quality here
    DATAdir = os.path.join(NMLDB_SIMU_PATH, simulation.indexingPath)
    print('Analyzing: ', DATAdir)

    #Order Parameters 
    system_quality = {}
    for lipid1 in simulation.getLipids():
        print('')
        print('Evaluating order parameter quality of simulation data in ' + simulation.indexingPath)
        
        OP_data_lipid = {}
        #convert elements to float because in some files the elements are strings
        try:
            for key, value in simulation.OPdata[lipid1].items():
                OP_array = [float(x) for x in simulation.OPdata[lipid1][key][0]]  
                OP_data_lipid[key] = OP_array
        except:
            continue
        
        #go through file paths in simulation.readme['EXPERIMENT']
        fragment_qual_dict = {}
        data_dict = {}
        
        for doi, path in simulation.readme['EXPERIMENT']['ORDERPARAMETER'][lipid1].items():
            print('Evaluating '+ lipid1 + ' lipid using experimental data from ' + doi + ' in ../../Data/experiments/OrderParameters/' + path)
                
            #load mapping file
            mapping_file = simulation.readme['COMPOSITION'][lipid1]['MAPPING']
            
            print(doi)
            OP_qual_data = {}
            # get readme file of the experiment
            experimentFilepath = os.path.join(NMLDB_EXP_PATH, "OrderParameters", path)
            print('Experimental data available at ' + experimentFilepath)
                
            READMEfilepathExperiment  = os.path.join(experimentFilepath, 'README.yaml')
            experiment = qq.Experiment()
            with open(READMEfilepathExperiment) as yaml_file_exp:
                readmeExp = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)
                experiment.readme = readmeExp

            exp_OP_filepath = os.path.join(experimentFilepath, lipid1 + '_Order_Parameters.json')
            lipidExpOPdata = {}
            try:
                with open(exp_OP_filepath) as json_file:
                    lipidExpOPdata = json.load(json_file)
            except FileNotFoundError:
                print("Experimental order parameter data do not exist for lipid " + lipid1 + ".")
                continue

            exp_error = 0.02

            for key in OP_data_lipid.keys():
                OP_array = OP_data_lipid[key].copy()
                try:
                    OP_exp = lipidExpOPdata[key][0][0]
                except KeyError:
                    continue
                else:
                    if not np.isnan(OP_exp):
                        OP_sim = OP_array[0]
                        op_sim_STEM = OP_array[2]
                        #changing to use shitness(TM) scale. This code needs to be cleaned
                        op_quality = qq.prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_STEM)
                        OP_array.append(OP_exp)
                        OP_array.append(exp_error)   #hardcoded!!!! 0.02 for all experiments
                        OP_array.append(op_quality)
                OP_qual_data[key] = OP_array    
                
            # save qualities of simulation compared to an experiment into a dictionary
            data_dict[doi] = OP_qual_data
                
            # calculate quality for molecule fragments headgroup, sn-1, sn-2
            fragments = qq.getFragments(mapping_file)
            fragment_qual_dict[doi] = qq.fragmentQuality(fragments, lipidExpOPdata, OP_data_lipid)
                
        try:
            fragment_quality_output = qq.fragmentQualityAvg(lipid1,fragment_qual_dict,fragments)
        except:
            print('no fragment quality')
            fragment_quality_output = {}

        try:
            system_quality[lipid1] = fragment_quality_output
        except:
            print('no system quality')
            system_quality[lipid1] = {}

        fragment_quality_file = os.path.join(DATAdir, lipid1 + '_FragmentQuality.json')

        FGout = False
        for FG in fragment_quality_output:
            #print(FG,fragment_quality_output[FG])
            if np.isnan(fragment_quality_output[FG]):
                continue
            if fragment_quality_output[FG] > 0:
                FGout = True
        if FGout:
            with open(fragment_quality_file, 'w') as f:  #write fragment qualities into a file for a molecule
                 json.dump(fragment_quality_output,f)

        # write into the OrderParameters_quality.json quality data file   
        outfile1 = os.path.join(DATAdir, lipid1 + '_OrderParameters_quality.json')
        try:
            with open(outfile1, 'w') as f:
                json.dump(data_dict,f, cls=CompactJSONEncoder)
        except:
            pass

    system_qual_output = qq.systemQuality(system_quality, simulation)
    # make system quality file

    outfile2 = os.path.join(DATAdir, 'SYSTEM_quality.json')
    SQout = False
    for SQ in system_qual_output:
        if system_qual_output[SQ] > 0:
            SQout = True
    if SQout:
        with open(outfile2, 'w') as f:
            json.dump(system_qual_output,f)
        print('Order parameter quality evaluated for '  + simulation.indexingPath)
        EvaluatedOPs += 1
        print('')
        
###################################################################################################################                
    #Form factor quality
        
    expFFpath = simulation.readme['EXPERIMENT']['FORMFACTOR']
    expFFdata = {}
    if len(expFFpath) > 0:
        expFFpath_full = os.path.join(NMLDB_EXP_PATH, "FormFactors", expFFpath)
        for subdir, dirs, files in os.walk( expFFpath_full ):
            for filename in files:  
                filepath = os.path.join(expFFpath_full, filename)
                if filename.endswith('.json'):
                    with open(filepath) as json_file:
                        expFFdata = json.load(json_file)
    
    simFFdata = simulation.FFdata

    if len(expFFpath) > 0 and len(simFFdata) > 0:
        ffQuality = qq.formfactorQuality(simFFdata, expFFdata)
        outfile3 = os.path.join(DATAdir, 'FormFactorQuality.json')
        with open(outfile3,'w') as f:
            json.dump(ffQuality,f)
        EvaluatedFFs += 1
        print('Form factor quality evaluated for ', DATAdir)
    else:
        ffQuality = 0


print('The number of systems with evaluated order parameters:', EvaluatedOPs)
print('The number of systems with evaluated form factors:', EvaluatedFFs)