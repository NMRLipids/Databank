import os
import yaml
import json
import matplotlib.pyplot as plt
import numpy as np
import math

from matplotlib import cm
from scipy.stats import norm

from databankLibrary import lipids_dict

databank_path = '../../Data/Simulations'


lipid_numbers_list =  lipids_dict.keys() # should contain all lipid names
ions_list = ['POT', 'SOD', 'CLA', 'CAL'] # should contain names of all ions

class Data:
    def __init__(self, molecule, data_path):
        self.molecule = molecule
        self.data = {}
        self.__load_data__(data_path)
    
    def __load_data__(self,data_path):
        with open(data_path) as json_file:
            self.data = json.load(json_file)
        json_file.close()

        
class Simulation:
    def __init__(self, readme, data, indexingPath):
        self.readme = readme
        self.data = data #object Data
        self.indexingPath = indexingPath
    
    def getLipids(self, molecules=lipid_numbers_list):
        lipids = []
        
        for key in self.readme['COMPOSITION'].keys():
            if key in molecules:
                lipids.append(key)

        return lipids
    
    def getIons(self, ions):
        simIons = []
        for key in self.readme['COMPOSITION'].keys():
            if key in ions:
                simIons.append(key)
                
        return simIons

    #fraction of each lipid with respect to total amount of lipids      
    def molarFraction(self, molecule,molecules=lipid_numbers_list): #only for lipids
        sum_lipids = 0
        number = sum(self.readme['COMPOSITION'][molecule]['COUNT']) 
        
        for key in self.readme['COMPOSITION'].keys():
            if key in molecules:
                sum_lipids += sum(self.readme['COMPOSITION'][key]['COUNT'])

        return number / sum_lipids
        
        
    # concentration of other molecules than lipids 
    # change name to ionConcentration()
    def ionConcentration(self, molecule, exp_counter_ions):
        lipids1 = self.getLipids()
        c_water = 55.5
        N_water = self.readme['COMPOSITION']['SOL']['COUNT']
        try:
            N_molecule = self.readme['COMPOSITION'][molecule]['COUNT'] #number of ions
        except KeyError:
            N_molecule = 0
        

        lipids2 = []
        if exp_counter_ions and N_molecule != 0:
            for lipid in lipids1:
                if molecule in exp_counter_ions.keys() and lipid == exp_counter_ions[molecule]:
                    N_lipid = self.readme['COMPOSITION'][lipid]['COUNT']
                 #   print(molecule + " " + lipid)
                 #   print(self.readme)
                    lipids2.append(sum(N_lipid))
        
        N_molecule = N_molecule - sum(lipids2)
       # print(N_molecule)
        
        c_molecule = (N_molecule * c_water) / N_water
        #print(c_molecule)
        
        return c_molecule
        
    def totalLipidConcentration(self):
        lipids = self.getLipids()
        c_water = 55.5
        N_water = self.readme['COMPOSITION']['SOL']['COUNT']
        N_lipids = 0
        for lipid in lipids:
            N_lipids += sum(self.readme['COMPOSITION'][lipid]['COUNT'])
        try:
            if (N_water / N_lipids) > 25 :
                tot_lipid_c = 'full hydration'
              #  print('full hydration')
            else:
                tot_lipid_c = (N_lipids * c_water) / N_water
        except ZeroDivisionError:
            print(self.readme)    
        return tot_lipid_c
        
##################
class Experiment:
    def __init__(self, readme, data, dataPath,exptype):
        self.readme = readme
        self.data = data #object Data
        self.dataPath = dataPath
        self.exptype = exptype
        
    def getLipids(self, molecules=lipid_numbers_list):
        lipids = []
        for key in molecules:
            try:
                if key in self.readme['MOLAR_FRACTIONS'].keys():
                    lipids.append(key)
            except KeyError:
                continue
        return lipids
    
    def getIons(self, ions):
        expIons = []
        
        for key in ions:
            try:
                if self.readme['ION_CONCENTRATIONS'][key] != 0: 
                    expIons.append(key)
            except KeyError:
                continue    
            try:
                if key in self.readme['COUNTER_IONS'].keys():
                    expIons.append(key)
            except AttributeError:
                continue
        return expIons
        
def loadSimulations():
    simulations = []
    for subdir, dirs, files in os.walk(r'../../Data/Simulations/'): 
        for filename1 in files:
            filepath = subdir + os.sep + filename1
        
            if filepath.endswith("README.yaml"):
                READMEfilepathSimulation = subdir + '/README.yaml'
                with open(READMEfilepathSimulation) as yaml_file_sim:
                    readmeSim = yaml.load(yaml_file_sim, Loader=yaml.FullLoader)
                    indexingPath = "/".join(filepath.split("/")[4:8])
                    simOPdata = [] #order parameter files for each type of lipid
                    for filename2 in files:
                        if filename2.endswith('OrderParameters.json'):
                            lipid_name = filename2.replace('OrderParameters.json', '')
                            dataPath = subdir + '/' + filename2
                            op_data = Data(lipid_name, dataPath)
                            simOPdata.append(op_data)
                    simulations.append(Simulation(readmeSim, simOPdata, indexingPath))
                    yaml_file_sim.close()    
    return simulations
 
def loadExperiments(experimentType):
#loop over the experiment entries in the experiment databank and read experiment readme and order parameter files into objects 
    experiments = []
    dataFile = ""
    if experimentType == 'OrderParameters':
        dataFile = 'Order_Parameters.json'
    elif experimentType == 'FormFactors':
        dataFile = 'form_factors.fourierFromFinalDensity.json'
    
    path = r'../../Data/experiments/'+experimentType+'/'
    for exp_subdir, exp_dirs, exp_files in os.walk(path):
        for filename1 in exp_files:
            filepath_exp = exp_subdir + os.sep + filename1
            if filepath_exp.endswith("README.yaml"):
                READMEfilepathExperiment = exp_subdir + '/README.yaml'
                #print(READMEfilepathExperiment)
                with open(READMEfilepathExperiment) as yaml_file_exp:
                    readmeExp = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)
                    for filename2 in exp_files:
                        if filename2.endswith(dataFile):
                            molecule_name = ""
                            dataPath = exp_subdir + '/' + filename2 #json!
                            if dataFile == "OrderParameters":
                                molecule_name = filename2.replace(dataFile,'')
                            elif dataFile == "FormFactors": 
                                molecule_name = 'system'   
                            #print(lipid_name)    
                            expOPdata = Data(molecule_name, dataPath)
                            experiments.append(Experiment(readmeExp, expOPdata, exp_subdir,experimentType))
                yaml_file_exp.close()
           
    return experiments
    
def findPairs(experiments):
    pairs = []
    for simulation in simulations:
        sim_lipids = simulation.getLipids()
        sim_total_lipid_concentration = simulation.totalLipidConcentration() 
        sim_ions = simulation.getIons(ions_list)
        t_sim = simulation.readme['TEMPERATURE']
    
        #calculate molar fractions from simulation
        sim_molar_fractions = {}
        for lipid in sim_lipids:
            sim_molar_fractions[lipid] = simulation.molarFraction(lipid)
        
        
        for experiment in experiments: 
            # print(experiment.readme)
            # check lipid composition matches the simulation
            exp_lipids = experiment.getLipids() 
            #  print(experiment.getLipids())

            exp_total_lipid_concentration = experiment.readme['TOTAL_LIPID_CONCENTRATION']
            exp_ions = experiment.getIons(ions_list)
            # print('experiment ions: ' + str(exp_ions)) 
            exp_counter_ions = experiment.readme['COUNTER_IONS']
        
            # calculate simulation ion concentrations
            sim_concentrations = {}
            for molecule in ions_list:
                sim_concentrations[molecule] = simulation.ionConcentration(molecule, exp_counter_ions)

        
            # continue if lipid compositions are the same
            if set(sim_lipids) == set(exp_lipids):
                # compare molar fractions
                mf_ok = 0
                for key in sim_lipids:
                    if (experiment.readme['MOLAR_FRACTIONS'][key] >= sim_molar_fractions[key] - 0.05) and (experiment.readme['MOLAR_FRACTIONS'][key] <= sim_molar_fractions[key]+ 0.05):
                        mf_ok +=1 

                c_ok = 0

                # compare ion concentrations 
                if set(sim_ions) == set(exp_ions):
                  #  print(sim_ions)
                  #  print(exp_ions)
                    for key in sim_ions:
                        if (experiment.readme['ION_CONCENTRATIONS'][key] >= sim_concentrations[key] - 0.05) and (experiment.readme['ION_CONCENTRATIONS'][key] <= sim_concentrations[key] + 0.05): 
                            c_ok += 1 


                switch = 0
            
                if (type(exp_total_lipid_concentration) == float) and (type(sim_total_lipid_concentration) == float): 
                    if ((exp_total_lipid_concentration >= sim_total_lipid_concentration - 0.1) and (exp_total_lipid_concentration <= sim_total_lipid_concentration + 0.1)):
                        switch = 1
                elif type(exp_total_lipid_concentration) == str and type(sim_total_lipid_concentration) == str:
                    if exp_total_lipid_concentration == sim_total_lipid_concentration:
                        switch = 1
                        
                if switch == 1:
            
                    #check temperature +/- 2 degrees
                    t_exp = experiment.readme['TEMPERATURE']
                    if (mf_ok == len(sim_lipids)) and (c_ok == len(sim_ions)) and (t_exp >= float(t_sim) - 2.0) and (t_exp <= float(t_sim) + 2.0):
                        #  print(simulation.indexingPath)
                        pairs.append([simulation, experiment])
                        # print(simulation.readme['SYSTEM'])
                        print(simulation.indexingPath)
                        print(experiment.dataPath)
                    
                    #Add path to experiment into simulation README.yaml
                    #many experiment entries can match to same simulation
                        exp_doi = experiment.readme['DOI'] 
                        exp_path = "/".join(experiment.dataPath.split("/")[5:8])
                        if experiment.exptype == "OrderParameters":
                            lipid = experiment.data.molecule
                            simulation.readme['EXPERIMENT']['ORDERPARAMETER'][lipid][exp_doi] = exp_path
                            print(simulation.readme['DOI'])
                            print(simulation.readme['EXPERIMENT'])
                            print('\n')
                        elif experiment.exptype == "FormFactors":
                            simulation.readme['EXPERIMENT']['FORMFACTOR'][exp_doi]=exp_path
                    
        outfileDICT = '../../Data/Simulations/'+ simulation.indexingPath + '/README.yaml'
    
        with open(outfileDICT, 'w') as f:
            yaml.dump(simulation.readme,f, sort_keys=False)
        f.close()
        
        return pairs
                        
##############################################
#loop over the simulations in the simulation databank and read simulation readme and order parameter files into objects
simulations = loadSimulations()
for simulation in simulations:
    simulation.readme['EXPERIMENT'] = {}
    for lipid in simulation.getLipids():
        simulation.readme['EXPERIMENT'][lipid] = {}
        
    outfileDICT = '../../Data/Simulations/'+ simulation.indexingPath + '/README.yaml'
    with open(outfileDICT, 'w') as f:
        yaml.dump(simulation.readme,f, sort_keys=False)

#load experiments
experimentsOrderParameters = loadExperiments('OrderParameters')
#experimentsFormFactors = loadExperiments('FormFactors')

# make empty dictionary for saving paths to matching experimental data entries                

                

                   

#Pair each simulation with an experiment with the closest matching temperature and composition
pairs = findPairs(experimentsOrderParameters)


print(len(experimentsOrderParameters))
print(len(simulations))

print("Found " + str(len(pairs)) + " pairs")  
#for pair in pairs:
#    print(pair[0].readme)
#    print(pair[1].readme)                    
