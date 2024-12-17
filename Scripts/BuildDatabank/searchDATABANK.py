#!/usr/bin/env python3
# coding: utf-8

import os
import yaml

from tqdm import tqdm
from typing import List, IO

from DatabankLib import NMLDB_SIMU_PATH, NMLDB_EXP_PATH
from DatabankLib.core import initialize_databank
from DatabankLib.databankLibrary import lipids_dict

import logging
logger = logging.getLogger("__name__")

lipid_numbers_list = lipids_dict.keys()   # should contain all lipid names
ions_list = ['POT', 'SOD', 'CLA', 'CAL']  # should contain names of all ions

LIP_CONC_REL_THRESHOLD = 0.15   # relative acceptable error for determination
# of the hydration in ssNMR # noqa


class Simulation:

    readme: dict
    indexingPath: str

    def __init__(self, readme):
        self.readme = readme
        self.indexingPath = readme['path']

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

    # fraction of each lipid with respect to total amount of lipids (only for lipids!)
    def molarFraction(self, molecule, molecules=lipid_numbers_list) -> float:
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
            N_molecule = self.readme['COMPOSITION'][molecule]['COUNT']  # number of ions
        except KeyError:
            N_molecule = 0

        lipids2 = []
        if exp_counter_ions and N_molecule != 0:
            for lipid in lipids1:
                if (molecule in exp_counter_ions.keys() and
                        lipid == exp_counter_ions[molecule]):
                    N_lipid = self.readme['COMPOSITION'][lipid]['COUNT']
                    lipids2.append(sum(N_lipid))

        N_molecule = N_molecule - sum(lipids2)
        c_molecule = (N_molecule * c_water) / N_water

        return c_molecule

    def totalLipidConcentration(self):
        lipids = self.getLipids()
        c_water = 55.5
        N_water = self.readme['COMPOSITION']['SOL']['COUNT']
        N_lipids = 0
        for lipid in lipids:
            N_lipids += sum(self.readme['COMPOSITION'][lipid]['COUNT'])
        try:
            if (N_water / N_lipids) > 25:
                tot_lipid_c = 'full hydration'
            else:
                tot_lipid_c = (N_lipids * c_water) / N_water
        except ZeroDivisionError:
            logger.warning("Division by zero when determining lipid concentration!")
            print(self.readme)
        return tot_lipid_c

##################


class Experiment:
    def __init__(self, readme: dict, molname: str, dataPath: str, exptype: str):
        self.readme = readme
        self.molname = molname    # <- the dictionary about existence of data files
        self.dataPath = dataPath  # .... for particular lipids
        self.exptype = exptype

    def getLipids(self, molecules=lipid_numbers_list) -> List[str]:
        lipids: List[str] = []
        for key in molecules:
            try:
                if key in self.readme['MOLAR_FRACTIONS'].keys():
                    lipids.append(key)
            except KeyError:
                continue
        return lipids

    def getIons(self, ions) -> List[str]:
        expIons: List[str] = []

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


def loadSimulations() -> List[Simulation]:
    """
    Generates the list of Simulation objects. Go through all README.yaml files.
    """

    systems = initialize_databank()
    simulations: List[Simulation] = []

    for system in systems:
        # conditions of exclusions
        try:
            if system['WARNINGS']['NOWATER']:
                continue
        except (KeyError, TypeError):
            pass

        simulations.append(Simulation(system))

    return simulations


def loadExperiments(experimentType: str) -> List[Experiment]:
    """
    Loops over the experiment entries in the experiment databank and read experiment
    readme and order parameter files into objects.
    """

    if experimentType == 'OrderParameters':
        dataFile = '_Order_Parameters.json'
    elif experimentType == 'FormFactors':
        dataFile = '_FormFactor.json'
    else:
        raise NotImplementedError(
            "Only OrderParameters and FormFactors types are implemented.")

    print("Build experiments [%s] index..." % experimentType, end='')
    rmIdx = []

    path = os.path.join(NMLDB_EXP_PATH, experimentType)
    for subdir, dirs, files in os.walk(path):
        for fn in files:
            if fn == 'README.yaml':
                rmIdx.append(subdir)
    print('%d READMEs loaded.' % len(rmIdx))

    print("Loading data for each experiment.")
    experiments: List[Experiment] = []
    for subdir in tqdm(rmIdx, desc='Experiment'):
        READMEfilepathExperiment = os.path.join(subdir, 'README.yaml')
        with open(READMEfilepathExperiment) as yaml_file_exp:
            readmeExp = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)

        for fname in os.listdir(subdir):
            if fname.endswith(dataFile):
                molecule_name = ""
                if experimentType == "OrderParameters":
                    molecule_name = fname.replace(dataFile, '')
                elif experimentType == "FormFactors":
                    molecule_name = 'system'
                experiments.append(Experiment(
                    readmeExp, molecule_name, subdir, experimentType))

    return experiments


def findPairs(experiments: List[Experiment], simulations: List[Simulation]):
    pairs = []
    for simulation in tqdm(simulations, desc='Simulation'):
        sim_lipids = simulation.getLipids()
        sim_total_lipid_concentration = simulation.totalLipidConcentration()
        sim_ions = simulation.getIons(ions_list)
        t_sim = simulation.readme['TEMPERATURE']

        # calculate molar fractions from simulation
        sim_molar_fractions = {}
        for lipid in sim_lipids:
            sim_molar_fractions[lipid] = simulation.molarFraction(lipid)

        for experiment in experiments:

            # check lipid composition matches the simulation
            exp_lipids = experiment.getLipids()

            exp_total_lipid_concentration = \
                experiment.readme['TOTAL_LIPID_CONCENTRATION']
            exp_ions = experiment.getIons(ions_list)
            exp_counter_ions = experiment.readme['COUNTER_IONS']

            # calculate simulation ion concentrations
            sim_concentrations = {}
            for molecule in ions_list:
                sim_concentrations[molecule] = simulation.ionConcentration(
                    molecule, exp_counter_ions)

            # continue if lipid compositions are the same
            if set(sim_lipids) == set(exp_lipids):
                # compare molar fractions
                mf_ok = 0
                for key in sim_lipids:
                    if ((experiment.readme['MOLAR_FRACTIONS'][key] >=
                         sim_molar_fractions[key] - 0.03) and
                        (experiment.readme['MOLAR_FRACTIONS'][key] <=
                         sim_molar_fractions[key] + 0.03)):
                        mf_ok += 1

                # compare ion concentrations
                c_ok = 0
                if set(sim_ions) == set(exp_ions):
                    for key in sim_ions:
                        if ((experiment.readme['ION_CONCENTRATIONS'][key] >=
                             sim_concentrations[key] - 0.05) and
                            (experiment.readme['ION_CONCENTRATIONS'][key] <=
                             sim_concentrations[key] + 0.05)):
                            c_ok += 1

                switch = 0

                if (isinstance(exp_total_lipid_concentration, (int, float)) and
                        isinstance(sim_total_lipid_concentration, (int, float))):
                    if ((exp_total_lipid_concentration / sim_total_lipid_concentration >
                         1 - LIP_CONC_REL_THRESHOLD) and
                        (exp_total_lipid_concentration / sim_total_lipid_concentration <
                         1 + LIP_CONC_REL_THRESHOLD)):
                        switch = 1
                elif ((type(exp_total_lipid_concentration) is str) and
                      (type(sim_total_lipid_concentration) is str)):
                    if exp_total_lipid_concentration == sim_total_lipid_concentration:
                        switch = 1

                if switch:
                    # check temperature +/- 2 degrees
                    t_exp = experiment.readme['TEMPERATURE']

                    if ((mf_ok == len(sim_lipids)) and
                        (c_ok == len(sim_ions)) and
                        (t_exp > float(t_sim) - 2.5) and
                            (t_exp < float(t_sim) + 2.5)):
                        # !we found the match!
                        pairs.append([simulation, experiment])

                        # Add path to experiment into simulation README.yaml
                        # many experiment entries can match to same simulation
                        exp_doi = experiment.readme['DOI']
                        exp_path = os.path.relpath(
                            experiment.dataPath,
                            start=os.path.join(NMLDB_EXP_PATH, experiment.exptype))
                        if experiment.exptype == "OrderParameters":
                            lipid = experiment.molname
                            simulation.readme['EXPERIMENT']['ORDERPARAMETER'][lipid][exp_doi] = exp_path # noqa
                        elif experiment.exptype == "FormFactors":
                            simulation.readme['EXPERIMENT']['FORMFACTOR'] = exp_path
                    else:
                        continue

        # sorting experiment lists to keep experimental order strict
        for _lipid in simulation.readme['EXPERIMENT']['ORDERPARAMETER'].keys():
            unsortDict = simulation.readme['EXPERIMENT']['ORDERPARAMETER'][_lipid].copy(
            )
            if not len(unsortDict):
                continue
            sortDict = dict(sorted(unsortDict.items()))
            simulation.readme['EXPERIMENT']['ORDERPARAMETER'][_lipid] = sortDict.copy()

        outfileDICT = os.path.join(
            NMLDB_SIMU_PATH, simulation.indexingPath, 'README.yaml')
        with open(outfileDICT, 'w') as f:
            if "path" in simulation.readme.keys():
                del (simulation.readme['path'])
            yaml.dump(simulation.readme, f, sort_keys=False, allow_unicode=True)

    return pairs


def logPairs(pairs, fd: IO[str]) -> None:
    """
    Write found correspondences into log file.

    pairs: [(Simulation, Experiment), ...]
    fd: file descriptor for writting into
    """

    for p in pairs:
        sim: Simulation = p[0]
        exp: Experiment = p[1]

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
    # TODO: check if EXPERIMENT section changed and trigger the action!
    for simulation in simulations:
        simulation.readme['EXPERIMENT'] = {}
        simulation.readme['EXPERIMENT']['ORDERPARAMETER'] = {}
        simulation.readme['EXPERIMENT']['FORMFACTOR'] = {}
        for lipid in simulation.getLipids():
            simulation.readme['EXPERIMENT']['ORDERPARAMETER'][lipid] = {}

        outfileDICT = os.path.join(
            NMLDB_SIMU_PATH, simulation.indexingPath, 'README.yaml')
        with open(outfileDICT, 'w') as f:
            yaml.dump(simulation.readme, f, sort_keys=False, allow_unicode=True)

    experimentsOrderParameters = loadExperiments('OrderParameters')
    experimentsFormFactors = loadExperiments('FormFactors')

    # Pair each simulation with an experiment with the closest matching temperature
    # and composition
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
