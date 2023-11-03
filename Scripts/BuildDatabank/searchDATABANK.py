import os
import yaml
import json

from databankLibrary import lipids_dict
from tqdm import tqdm

databank_path = '../../Data/Simulations'
expbank_path = '../../Data/experiments'

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

        
class Simulation:
    def __init__(self, readme, data, indexingPath):
        self.readme = readme
        self.data = {}      
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
    """
    Generates the list of Simulation objects. Go through all README.yaml files.
    """

    print("Build simulation tree index...", end='')
    rmIdx = []
    for subdir, dirs, files in os.walk(databank_path): 
        for fn in files:
            if fn == 'README.yaml':
                rmIdx.append(subdir)
    print('%d READMEs loaded.' % len(rmIdx) )

    print("Loading information about every simulation in the bank.")
    simulations = []
    for subdir in tqdm(rmIdx, "Simulations"):
        READMEfilepathSimulation = os.path.join(subdir, 'README.yaml')
        # it exists because we collected only dirs with README.yaml
        with open(READMEfilepathSimulation) as yaml_file_sim:
            readmeSim = yaml.load(yaml_file_sim, Loader=yaml.FullLoader)
        try:
            if readmeSim['WARNINGS']['NOWATER']:
                continue
        except:
            pass
        indexingPath = os.path.relpath(subdir, start=databank_path)
        simOPdata = [] # order parameter files for each type of lipid
        simData = {}

        for filename in os.listdir(subdir): 
            filepath = os.path.join(subdir, filename)
        
            if filename == 'fourierFromFinalDensity.json': ## TODO: WHAT IS THAT???
                simData['FormFactors'] = Data("system", filepath)
            elif filename.endswith('OrderParameters.json'):
                lipid_name = filename.replace('OrderParameters.json', '')
            #  print(dataPath)
                op_data = Data(lipid_name, filepath)
                simOPdata.append(op_data)
                            
        simData['OrderParameters'] = simOPdata
        simulations.append(Simulation(readmeSim, simData, indexingPath))

    return simulations
 
def loadExperiments(experimentType):
    """
    Loops over the experiment entries in the experiment databank and read experiment 
    readme and order parameter files into objects.
    """

    if experimentType == 'OrderParameters':
        dataFile = '_Order_Parameters.json'
    elif experimentType == 'FormFactors':
        dataFile = '_FormFactor.json'
    else:
        raise NotImplementedError("Only OrderParameters and FormFactors types are implemented.")

    print("Build experiments [%s] index..." % experimentType, end='')
    rmIdx = []

    path = os.path.join(expbank_path, experimentType)
    for subdir, dirs, files in os.walk(path): 
        for fn in files:
            if fn == 'README.yaml':
                rmIdx.append(subdir)
    print('%d READMEs loaded.' % len(rmIdx) )

    print("Loading data for each experiment.")
    experiments = []
    for subdir in tqdm(rmIdx, desc='Experiment'):
        READMEfilepathExperiment = os.path.join(subdir, 'README.yaml')
        with open(READMEfilepathExperiment) as yaml_file_exp:
            readmeExp = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)
        opData = {}

        for fname in os.listdir(subdir):
            dataPath = os.path.join(subdir, fname)
            if fname.endswith(dataFile):
                molecule_name = ""
                if experimentType == "OrderParameters":
                    molecule_name = fname.replace(dataFile,'')
                elif experimentType == "FormFactors": 
                    molecule_name = 'system'   
                expData = Data(molecule_name, dataPath)
                experiments.append(Experiment(readmeExp, expData, subdir, experimentType))
          
    return experiments
 
def findPairs(experiments, simulations):
    pairs = []
    for simulation in tqdm(simulations, desc='Simulation'):
        sim_lipids = simulation.getLipids()
        sim_total_lipid_concentration = simulation.totalLipidConcentration() 
        sim_ions = simulation.getIons(ions_list)
        t_sim = simulation.readme['TEMPERATURE']
        
        #calculate molar fractions from simulation
        sim_molar_fractions = {}
        for lipid in sim_lipids:
            sim_molar_fractions[lipid] = simulation.molarFraction(lipid)

        for experiment in experiments: 

            # check lipid composition matches the simulation
            exp_lipids = experiment.getLipids() 
            
            exp_total_lipid_concentration = experiment.readme['TOTAL_LIPID_CONCENTRATION']
            exp_ions = experiment.getIons(ions_list)
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
                    if ( (experiment.readme['MOLAR_FRACTIONS'][key] >= sim_molar_fractions[key] - 0.03) and 
                         (experiment.readme['MOLAR_FRACTIONS'][key] <= sim_molar_fractions[key]+ 0.03) ):
                        mf_ok +=1 

                # compare ion concentrations 
                c_ok = 0
                if set(sim_ions) == set(exp_ions):
                    for key in sim_ions:
                        if ( (experiment.readme['ION_CONCENTRATIONS'][key] >= sim_concentrations[key] - 0.05) and 
                             (experiment.readme['ION_CONCENTRATIONS'][key] <= sim_concentrations[key] + 0.05) ): 
                            c_ok += 1 

                switch = 0
            
                if ( (type(exp_total_lipid_concentration) == float) and 
                     (type(sim_total_lipid_concentration) == float) ):
                    if ( (exp_total_lipid_concentration >= sim_total_lipid_concentration - 0.1) and 
                         (exp_total_lipid_concentration <= sim_total_lipid_concentration + 0.1) ):
                        switch = 1
                elif ( (type(exp_total_lipid_concentration) == str) and 
                       (type(sim_total_lipid_concentration) == str) ):
                    if exp_total_lipid_concentration == sim_total_lipid_concentration:
                        switch = 1
                        
                if switch:
                    #check temperature +/- 2 degrees
                    t_exp = experiment.readme['TEMPERATURE']
                    
                    if ( (mf_ok == len(sim_lipids)) and 
                         (c_ok == len(sim_ions)) and 
                         (t_exp >= float(t_sim) - 2.0) and 
                         (t_exp <= float(t_sim) + 2.0) ):
                        # !we found the match!
                        pairs.append([simulation, experiment])

                        # Add path to experiment into simulation README.yaml
                        # many experiment entries can match to same simulation
                        exp_doi = experiment.readme['DOI'] 
                        exp_path = os.path.relpath(experiment.dataPath, start=os.path.join(expbank_path, experiment.exptype))
                        if experiment.exptype == "OrderParameters":
                            lipid = experiment.data.molecule
                            simulation.readme['EXPERIMENT']['ORDERPARAMETER'][lipid][exp_doi] = exp_path
                        elif experiment.exptype == "FormFactors":
                            simulation.readme['EXPERIMENT']['FORMFACTOR'] = exp_path
                    else:
                        continue
        
        # sorting experiment lists to keep experimental order strict
        for _lipid in simulation.readme['EXPERIMENT']['ORDERPARAMETER'].keys():
            unsortDict = simulation.readme['EXPERIMENT']['ORDERPARAMETER'][_lipid].copy()
            if not len(unsortDict):
                continue
            sortDict = dict(sorted(unsortDict.items()))
            simulation.readme['EXPERIMENT']['ORDERPARAMETER'][_lipid] = sortDict.copy()

        outfileDICT = os.path.join(databank_path, simulation.indexingPath, 'README.yaml')
        with open(outfileDICT, 'w') as f:
            yaml.dump(simulation.readme, f, sort_keys=False)
        
    return pairs

def logPairs(pairs, fd):
    """
    Write found correspondences into log file.

    pairs: [(Simulation, Experiment), ...]
    fd: file descriptor for writting into
    """

    for p in pairs:
        sim : Simulation = p[0]
        exp : Experiment = p[1]

        sysn = sim.readme['SYSTEM']
        simp = sim.indexingPath

        expp = exp.dataPath
        expd = exp.readme['DOI']

        fd.write(f"""
--------------------------------------------------------------------------------
Simulation:
 - {sysn}
 - {simp}
Experiment:
 - {expd}
 - {expp}""")
        # end for
    fd.write("""
--------------------------------------------------------------------------------
    \n""")

def main():
    """
    Main program function. Not for exporting.
    """

    simulations = loadSimulations()

    # clear all EXPERIMENT sections in all simulations
    for simulation in simulations:
        simulation.readme['EXPERIMENT'] = {}
        simulation.readme['EXPERIMENT']['ORDERPARAMETER']= {}
        simulation.readme['EXPERIMENT']['FORMFACTOR']= {}
        for lipid in simulation.getLipids():
            simulation.readme['EXPERIMENT']['ORDERPARAMETER'][lipid] = {}
        
        outfileDICT = os.path.join(databank_path, simulation.indexingPath, 'README.yaml')
        with open(outfileDICT, 'w') as f:
            yaml.dump(simulation.readme, f, sort_keys=False)

    experimentsOrderParameters = loadExperiments('OrderParameters')
    experimentsFormFactors = loadExperiments('FormFactors')

    # Pair each simulation with an experiment with the closest matching temperature and composition
    with open('search-databank-pairs.log', 'w') as logf:
        print("Scanning simulation-experiment pairs among order parameter experiments.")
        pairsOP = findPairs(experimentsOrderParameters, simulations)
        logf.write("=== OP PAIRS ===\n")
        logPairs(pairsOP, logf)
        print("Scanning simulation-experiment pairs among form factor experiments.")
        pairsFF = findPairs(experimentsFormFactors, simulations)
        logf.write("=== FF PAIRS ===\n")
        logPairs(pairsFF, logf)

    '''
    for pair in pairsFF:
        print('#################')
        print(pair[0].readme)
        print(pair[0].indexingPath)
        print("#")
        print(pair[1].readme)         
    '''

    print("Found order parameter data for " + str(len(pairsOP)) + " pairs")  
    print("Found form factor data for " + str(len(pairsFF)) + " pairs")

if __name__ == "__main__":
    main()