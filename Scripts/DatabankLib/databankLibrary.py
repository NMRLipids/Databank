# Library that contains all dictionaries and functions used in building and analyzing the NMR lipids databank

####################DICTIONARIES#####################################
# Define dictionaries
#
# Dictionary of lipids.
#
# If you add a lipid which is not yet in the databank, you have to add it here

import copy
import hashlib
import urllib, logging

import matplotlib.pyplot as plt
from tqdm import tqdm
from typing import List

import json
import sys


__pdoc__ = {}

logger = logging.getLogger("__name__")

from .databank_defs import *
from .databankio import resolve_download_file_url

##############CLASS FOR LOOPING OVER SYSTEMS#######################################

import yaml


class databank:
    """ :meta private: 
    Representation of all simulation in the NMR lipids databank. 

        `path` should be the local location of /Data/Simulations/ in the NMRlipids databank folder. Example usage to loop over systems: 
   
            path = '../../Data/Simulations/'
            db_data = databank(path)
            systems = db_data.get_systems()

            for system in systems:
                print(system)
  
    """
    
    def __init__(self, path=r"../../Data/Simulations/"):
        self.path = path
        self.systems = []
        self.__load_systems__(path)
        print('Databank initialized from the folder:', os.path.realpath(path))

    def __load_systems__(self, path):
        for subdir, dirs, files in os.walk(path):
            for filename in files:
                filepath = os.path.join(subdir, filename)
                #print(filepath)
                if filename == "README.yaml":
                    with open(filepath) as yaml_file:
                        content = yaml.load(yaml_file, Loader=yaml.FullLoader)
                        size = len(filepath)
                        sizePath = len(path)
                        content["path"] = filepath[sizePath : size - 11]
                        self.systems.append(content)

    def get_systems(self):
        """ Returns a list of all systems in the NMRlipids databank """
        return self.systems

#    def pie_temperature(self):
#        list_feature = [int(float(system["TEMPERATURE"])) for system in self.systems]
#        import collections

#        counter = collections.Counter(list_feature)
#        plt.pie(counter.values(), labels=counter.keys(), normalize=True)


#########################FUNCTIONS###################################################
# functions used in building and analyzing the databank

def initialize_databank(databankPath):
    """ 
    Intializes the NMRlipids databank.

    :param databankPath: path for the local location of the NMRlipids databank, for example ``../../Databank``
    :return: list of dictionaries that contain the content of README.yaml files for each system.  
    """
    #sys.path.insert(1, databankPath + '/Scripts/BuildDatabank/')
    #from databankLibrary import download_link, lipids_dict, databank
    path = databankPath + '/Data/Simulations/'
    db_data = databank(path)
    systems = db_data.get_systems()
    return systems


def print_README(system):
    """ 
    Prints the content of ``system`` dictionary in human readable format. 

    :param system: NMRlipids databank dictionary defining a simulation.

    """
    if system == 'example':
        readmePath = '../data/simulations/READMEexplanations.yaml'
        with open(readmePath, 'r') as file:
            readmeFile = yaml.safe_load(file)
    else:
        readmeFile = system
            
    for key in readmeFile:
        print('\033[1m' + key + ":" + '\033[0m')
        print(" ", readmeFile[key])



def CalcAreaPerMolecule(system):
    """ 
    Calculates average area per lipid for a simulation defined with ``system``. 
    It is using the ``apl.json`` file where area per lipid as a function of time calculated by the ``calcAPL.py`` is stored. 

    :param system: NMRlipids databank dictionary defining a simulation.

    :return: area per lipid (Å^2)
    """
    APLpath = os.path.dirname(os.path.realpath(__file__)) + '/../../Data/Simulations/' + system['path'] + 'apl.json'
    try:
        f = open(APLpath)
        APLdata = json.load(f)
        sumAPL = 0
        sumIND = 0
        for i,j in APLdata.items():
            sumAPL += j
            sumIND += 1
        APL = sumAPL/sumIND
        return(APL)
    except:
        print('apl.json not found from' + APLpath)
        

def GetThickness(system):
    """ 
    Gets thickness for a simulation defined with ``system`` from the ``thickness.json`` file 
    where thickness calculated by the ``calc_thickness.py`` is stored. 

    :param system: NMRlipids databank dictionary defining a simulation.

    :return: membrane thickess (nm)
    """
    ThicknessPath = os.path.dirname(os.path.realpath(__file__)) + '/../../Data/Simulations/' + system['path'] + 'thickness.json'
    try:
        f = open(ThicknessPath)
        thickness = json.load(f)
        return(thickness)
    except:
        return 0
        #print('thickness.json not found from' + system['path'])


def ShowEquilibrationTimes(system):
    """ 
    Prints relative equilibration time for each lipid within a simulation defined by ``system``. 
    Relative equilibration times are calculated with ``NMRPCA_timerelax.py`` and stored in ``eq_times.json`` files.

    :param system: NMRlipids databank dictionary defining a simulation.
    """
    
    EqTimesPath = os.path.dirname(os.path.realpath(__file__)) + '/../../Data/Simulations/' + system['path'] + 'eq_times.json'
    
    try:
        #if (os.path.isfile(EqTimesPath)):
        with open(EqTimesPath) as f:
            EqTimeDict = json.load(f)
    except:
        print('eq_times.json not found')
        exit()

    for i in EqTimeDict:
        print(i+':', EqTimeDict[i])
            
def GetNlipids(system):
    """ 
    Returns the total number of lipids in a simulation defined by ``system``. 

    :param system: NMRlipids databank dictionary defining a simulation.

    :return: the total number of lipids in the ``system``.
    """
    Nlipid = 0
    for molecule in system['COMPOSITION']:
        if molecule in lipids_dict:
            Nlipid += np.sum(system['COMPOSITION'][molecule]['COUNT'])
    return Nlipid

def getLipids(system, molecules=lipids_dict.keys()):
    """
    Returns a string using MDAnalysis notation that can used to select all lipids from the ``system``.

    :param system: NMRlipids databank dictionary defining a simulation.

    :return: a string using MDAnalysis notation that can used to select all lipids from the ``system``.
    """
    lipids = 'resname '

    for key in system['COMPOSITION'].keys():
        if key in molecules:
            m_file = system['COMPOSITION'][key]['MAPPING']
            with open('../../Databank/Scripts/BuildDatabank/mapping_files/'+m_file,"r") as yaml_file:
                mapping_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)
                #for line in f:
                #    if len(line.split()) > 2 and "Individual atoms" not in line:
            for atom in mapping_dict:
                print(mapping_dict[atom])
                try:
                    #        if line.split()[2] not in lipids:
                    lipids = lipids + mapping_dict[atom]['RESIDUE'] + ' or resname '
                except:
                    lipids = lipids + system['COMPOSITION'][key]['NAME'] + ' or resname '
                    break
                    #continue
                    #else:

                 
    lipids = lipids[:-12]
    print(lipids)
    return lipids

        
__pdoc__['plotFormFactor'] = False
def plotFormFactor(expFormFactor,k,legend,PlotColor):
    """:meta private:"""
    xValues = []
    yValues = []
    for i in expFormFactor:
        xValues.append(i[0])
        yValues.append(k*i[1])
    plt.plot(xValues,yValues,label = legend,color=PlotColor,linewidth=4.0)
    plt.xlabel(r'$q_{z} [Å^{-1}]$',size=20)
    plt.ylabel(r'$|F(q_{z})|$',size=20)
    plt.xticks(size=20)
    plt.yticks(size=20)
    #plt.yticks(color = 'w')
    plt.xlim([0,0.69])
    plt.ylim([-10,250])
    plt.legend(loc="upper right")
    plt.savefig('FormFactor.pdf')


__pdoc__['plotOrderParameters'] = False
def plotOrderParameters(OPsim, OPexp):
    """:meta private:"""
    xValuesHG = []
    xValuesSN1 = []
    xValuesSN2 = []
    
    yValuesHGsim = []
    yValuesSN1sim = []
    yValuesSN2sim = []
    yValuesHGsimERR = []
    yValuesSN1simERR = []
    yValuesSN2simERR = []
    yValuesHGexp = []
    yValuesSN1exp = []
    yValuesSN2exp = []
    xValuesHGexp = []
    xValuesSN1exp = []
    xValuesSN2exp = []

    sn1carbons = {'M_G1C3_M M_G1C3H1_M' :2,
                  'M_G1C3_M M_G1C3H2_M' :2,
                  'M_G1C4_M M_G1C4H1_M' :3,
                  'M_G1C4_M M_G1C4H2_M' :3,
                  'M_G1C5_M M_G1C5H1_M' :4,
                  'M_G1C5_M M_G1C5H2_M' :4,
                  'M_G1C6_M M_G1C6H1_M' :5,
                  'M_G1C6_M M_G1C6H2_M' :5,
                  'M_G1C7_M M_G1C7H1_M' :6,
                  'M_G1C7_M M_G1C7H2_M' :6,
                  'M_G1C8_M M_G1C8H1_M' :7,
                  'M_G1C8_M M_G1C8H2_M' :7,
                  'M_G1C9_M M_G1C9H1_M' :8,
                  'M_G1C9_M M_G1C9H2_M' :8,
                  'M_G1C10_M M_G1C10H1_M' :9,
                  'M_G1C10_M M_G1C10H2_M' :9,
                  'M_G1C11_M M_G1C11H1_M' :10,
                  'M_G1C11_M M_G1C11H2_M' :10,
                  'M_G1C12_M M_G1C12H1_M' :11,
                  'M_G1C12_M M_G1C12H2_M' :11,
                  'M_G1C13_M M_G1C13H1_M' :12,
                  'M_G1C13_M M_G1C13H2_M' :12,
                  'M_G1C14_M M_G1C14H1_M' :13,
                  'M_G1C14_M M_G1C14H2_M' :13,
                  'M_G1C15_M M_G1C15H1_M' :14,
                  'M_G1C15_M M_G1C15H2_M' :14,
                  'M_G1C16_M M_G1C16H1_M' :15,
                  'M_G1C16_M M_G1C16H2_M' :15,
                  'M_G1C17_M M_G1C17H1_M' :16,
                  'M_G1C17_M M_G1C17H2_M' :16,
                  'M_G1C17_M M_G1C17H3_M' :16,
                 }
    
    sn2carbons = {'M_G2C3_M M_G2C3H1_M' :2,
                  'M_G2C3_M M_G2C3H2_M' :2,
                  'M_G2C4_M M_G2C4H1_M' :3,
                  'M_G2C4_M M_G2C4H2_M' :3,
                  'M_G2C5_M M_G2C5H1_M' :4,
                  'M_G2C5_M M_G2C5H2_M' :4,
                  'M_G2C6_M M_G2C6H1_M' :5,
                  'M_G2C6_M M_G2C6H2_M' :5,
                  'M_G2C7_M M_G2C7H1_M' :6,
                  'M_G2C7_M M_G2C7H2_M' :6,
                  'M_G2C8_M M_G2C8H1_M' :7,
                  'M_G2C8_M M_G2C8H2_M' :7,
                  'M_G2C9_M M_G2C9H1_M' :8,
                  'M_G2C9_M M_G2C9H2_M' :8,
                  'M_G2C10_M M_G2C10H1_M' :9,
                  'M_G2C10_M M_G2C10H2_M' :9,
                  'M_G2C11_M M_G2C11H1_M' :10,
                  'M_G2C11_M M_G2C11H2_M' :10,
                  'M_G2C12_M M_G2C12H1_M' :11,
                  'M_G2C12_M M_G2C12H2_M' :11,
                  'M_G2C13_M M_G2C13H1_M' :12,
                  'M_G2C13_M M_G2C13H2_M' :12,
                  'M_G2C14_M M_G2C14H1_M' :13,
                  'M_G2C14_M M_G2C14H2_M' :13,
                  'M_G2C15_M M_G2C15H1_M' :14,
                  'M_G2C15_M M_G2C15H2_M' :14,
                  'M_G2C16_M M_G2C16H1_M' :15,
                  'M_G2C16_M M_G2C16H2_M' :15,
                  'M_G2C17_M M_G2C17H1_M' :16,
                  'M_G2C17_M M_G2C17H2_M' :16,
                  'M_G2C17_M M_G2C17H3_M' :16,
                  'M_G2C18_M M_G2C18H1_M' :17,
                  'M_G2C18_M M_G2C18H2_M' :17,
                  'M_G2C18_M M_G2C18H3_M' :17,
                  'M_G2C19_M M_G2C19H1_M' :18,
                  'M_G2C19_M M_G2C19H2_M' :18,
                  'M_G2C19_M M_G2C19H3_M' :18,
                 }
    
    HGcarbons = {'M_G3N6C1_M M_G3N6C1H1_M' : 1,
                 'M_G3N6C1_M M_G3N6C1H2_M' : 1,
                 'M_G3N6C1_M M_G3N6C1H3_M' : 1,
                 'M_G3N6C2_M M_G3N6C2H1_M' : 1,
                 'M_G3N6C2_M M_G3N6C2H2_M' : 1,
                 'M_G3N6C2_M M_G3N6C2H3_M' : 1,
                 'M_G3N6C3_M M_G3N6C3H1_M' : 1,
                 'M_G3N6C3_M M_G3N6C3H2_M' : 1,
                 'M_G3N6C3_M M_G3N6C3H3_M' : 1,
                 'M_G3C5_M M_G3C5H1_M' : 2,
                 'M_G3C5_M M_G3C5H2_M' : 2,
                 'M_G3C4_M M_G3C4H1_M' : 3,
                 'M_G3C4_M M_G3C4H2_M' : 3,
                 'M_G3_M M_G3H1_M' : 4,
                 'M_G3_M M_G3H2_M' : 4,
                 'M_G2_M M_G2H1_M' : 5,
                 'M_G1_M M_G1H1_M' : 6,
                 'M_G1_M M_G1H2_M' : 6,
                 }
    
    
    for key in OPsim:
        if 'M_G1C' in key:
            try:
                xValuesSN1.append(sn1carbons[key])
                yValuesSN1sim.append(float(OPsim[key][0][0]))
                yValuesSN1simERR.append(float(OPsim[key][0][2]))
                yValuesSN1exp.append(OPexp[key][0][0])
                xValuesSN1exp.append(sn1carbons[key])
            except:
                pass
        elif 'M_G2C' in key:
            try:
                xValuesSN2.append(sn2carbons[key])
                yValuesSN2sim.append(float(OPsim[key][0][0]))
                yValuesSN2simERR.append(float(OPsim[key][0][2]))
                yValuesSN2exp.append(OPexp[key][0][0])
                xValuesSN2exp.append(sn2carbons[key])
            except:
                pass
        elif 'M_G3' in key or 'M_G2_M' in key or 'M_G1_M' in key:
            try:
                xValuesHG.append(HGcarbons[key])
                yValuesHGsim.append(float(OPsim[key][0][0]))
                yValuesHGsimERR.append(float(OPsim[key][0][2]))
                yValuesHGexp.append(OPexp[key][0][0])
                xValuesHGexp.append(HGcarbons[key])
            except:
                pass
    #print(xValues,yValues)
    plt.rc('font', size=15)
    #plt.plot(xValuesHG,yValuesHGsim,'.',color='red',markersize=15)
    plt.errorbar(xValuesHGexp,yValuesHGexp, yerr = 0.02,fmt='.',color='black',markersize=25)
    plt.errorbar(xValuesHG,yValuesHGsim, yerr = yValuesHGsimERR,fmt='.',color='red',markersize=20)
    #plt.plot(xValuesHG,yValuesHGexp,'.',color='black',markersize=15)
    my_xticks = ['\u03B3','\u03B2','\u03B1','$g_{1}$','$g_{2}$','$g_{3}$']
    plt.xticks([1,2,3,4,5,6], my_xticks,size=20)
    #plt.xlabel('Carbon')
    #plt.ylim([-0.25,0.25])
    plt.yticks(size=20)
    #plt.yticks(color = 'w')
    plt.ylabel(r'$S_{CH}$',size=25)
    plt.savefig('HG.pdf')
    plt.show()
    
    plt.text(2, -0.04, 'sn-1', fontsize=25)
    plt.xticks(np.arange(min(xValuesSN1), max(xValuesSN1)+1, 2.0))
    plt.plot(xValuesSN1,yValuesSN1sim,color='red')
    plt.plot(xValuesSN1exp,yValuesSN1exp,color='black')
    plt.errorbar(xValuesSN1,yValuesSN1sim, yerr = yValuesSN1simERR,fmt='.',color='red',markersize=25)
    plt.errorbar(xValuesSN1exp,yValuesSN1exp, yerr = 0.02, fmt='.',color='black',markersize=20)
    #plt.xlabel('Carbon')
    plt.ylabel(r'$S_{CH}$',size=25)
    plt.xticks(size=20)
    plt.yticks(size=20)
    #plt.yticks(color = 'w')
    #plt.ylim([-0.3,0.01])
    plt.savefig('sn-1.pdf')
    plt.show()
    
    plt.text(2, -0.04, 'sn-2', fontsize=25)
    plt.xticks(np.arange(min(xValuesSN2), max(xValuesSN2)+1, 2.0))
    plt.plot(xValuesSN2,yValuesSN2sim,color='red')
    plt.plot(xValuesSN2exp,yValuesSN2exp,color='black')
    plt.errorbar(xValuesSN2,yValuesSN2sim,yValuesSN2simERR,fmt='.',color='red',markersize=25)
    plt.errorbar(xValuesSN2exp,yValuesSN2exp, yerr = 0.02, fmt='.',color='black',markersize=20)
    plt.xlabel('Carbon',size=25)
    plt.ylabel(r'$S_{CH}$',size=25)
    plt.xticks(size=20)
    plt.yticks(size=20)
    #plt.yticks(color = 'w')
    #plt.ylim([-0.25,0.01])
    plt.savefig('sn-2.pdf')
    plt.show()


def plotSimulation(ID, lipid):
    """
    Creates plots of form factor and C-H bond order parameters for the selected ``lipid`` from a simulation with the given ``ID`` number. 

    :param ID: NMRlipids databank ID number of the simulation
    :param lipid: universal molecul name of the lipid

    """
    DataBankPath =  os.path.dirname(os.path.realpath(__file__)) + '/../../'
    systems = initialize_databank(DataBankPath)
    #print(systems)
    for system in systems:
        #print(system)
        if system['ID'] == ID:
             path = DataBankPath + 'Data/Simulations/' + system['path']
    #lipid = 'POPC'
    FFpathSIM = path + 'FormFactor.json'
    OPpathSIM = path + lipid + 'OrderParameters.json'
    READMEfilepath = path + '/README.yaml'
    FFQualityFilePath = path + '/FormFactorQuality.json'

    
    with open(READMEfilepath) as yaml_file:
        readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
   
    print('DOI: ', readme['DOI'])
    
    try:
        with open(FFQualityFilePath) as json_file:
            FFq = json.load(json_file)
        print('Form factor quality: ', FFq[0])
        for subdir, dirs, files in os.walk(DataBankPath + 'Data/experiments/FormFactors/' + readme['EXPERIMENT']['FORMFACTOR'] + '/'):
            for filename in files:
                #filepath = '../../Data/experiments/FormFactors/' + expFFpath + '/' + filename
                if filename.endswith('_FormFactor.json'):
                    FFpathEXP = subdir + filename
        #FFpathEXP =  DataBankPath + 'experiments/FormFactors/' + readme['EXPERIMENT']['FORMFACTOR'] + '/POPS_ULV_25Cin0D_SHE_FormFactor.json'
        with open(FFpathEXP) as json_file:
            FFexp = json.load(json_file)
    
    except:
        print('Force field quality not found')
    
    
    with open(OPpathSIM) as json_file:
        OPsim = json.load(json_file)

    OPexp = {}
    for expOPfolder in list(readme['EXPERIMENT']['ORDERPARAMETER'][lipid].values()):
        #expOPfolder = list(readme['EXPERIMENT']['ORDERPARAMETER'][lipid].values())[0]
        OPpathEXP =  DataBankPath + 'Data/experiments/OrderParameters/' + expOPfolder + '/' + lipid + '_Order_Parameters.json'
        #print(OPpathEXP)
        with open(OPpathEXP) as json_file:
            OPexp.update(json.load(json_file))
    #print(OPexp)

    try:
        with open(FFpathSIM) as json_file:
            FFsim = json.load(json_file)
        plotFormFactor(FFsim,1, 'Simulation','red')
        plotFormFactor(FFexp,FFq[1], 'Experiment','black')
        plt.show()
    except:
        plt.show()
        print('Form factor plotting failed')
    
    plotOrderParameters(OPsim, OPexp)
    #print(OPsim)
    #print(OPexp)

        
def simulation2universal_atomnames(system,molecule,atom):
    """
    Get force field specific atom name corresponding to universal atom name from the ``system``.

    :param mapping_file: path for the mapping file
    :param atom1: universal atom name

    :return: force field specific atom name
    """

    mapping_file_path = os.path.dirname(os.path.realpath(__file__)) + '/mapping_files/' + system['COMPOSITION'][molecule]['MAPPING']
    with open(mapping_file_path, "rt") as mapping_file:
        mapping = yaml.load(mapping_file, Loader=yaml.FullLoader)
        try:
            m_atom1 = mapping[atom]["ATOMNAME"]
        except:
            print(atom, ' was not found from ', str(mapping_file_path))
            return   
    return m_atom1

def read_mapping_file(mapping_file, atom1):
    """
    Get force field specific atom name corresponding to universal atom name from the mapping file.

    :param mapping_file: path for the mapping file
    :param atom1: universal atom name

    :return: force field specific atom name
    """
    with open(mapping_file, "rt") as mapping_file:
        mapping = yaml.load(mapping_file, Loader=yaml.FullLoader)
        m_atom1 = mapping[atom1]["ATOMNAME"]
    return m_atom1

def loadMappingFile(mapping_file):
    """ 
    Load mapping file into a dictionary

    :param: name of the mapping file
    
    :return: mapping dictionary
    """
    mapping_dict = {}
    with open(os.path.dirname(os.path.realpath(__file__)) + '/mapping_files/' + mapping_file, "r") as yaml_file:
        mapping_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)
    yaml_file.close()
    
    return mapping_dict

def getAtoms(system, lipid):
    """
    Return system specific atom names of a lipid

    :param system: system dictionary
    :param lipid: universal lipid name

    :return: string of system specific atom names
    """
    
    atoms = ""
    path_to_mapping_file = system['COMPOSITION'][lipid]['MAPPING']
    mapping_dict = loadMappingFile(path_to_mapping_file)
    for key in mapping_dict:
        atoms = atoms + ' ' + mapping_dict[key]['ATOMNAME']
  
    return atoms

def getUniversalAtomName(system,atomName,lipid):
    """
    Returns the universal atom name corresponding the simulation specific ``atomName`` of a ``lipid`` in a simulation defined by the ``system``.

    :param system: system dictionary
    :param atomName: simulation specific atomname
    :param lipid: universal lipid name

    :return: universal atomname (string)
    """
    try:
        mappingFile = system['COMPOSITION'][lipid]['MAPPING']
    except:
        print('Mapping file was not found')
        return None

    mappingDict = loadMappingFile(mappingFile)

    for universalName in mappingDict:
        simName = mappingDict[universalName]['ATOMNAME']
        if simName == atomName:
            return universalName

    print('Atom was not found')
    return None


#######################ORDER PARAMETERS######################################
"""
    :meta private:
 calculation of order parameters of lipid bilayers
 from a MD trajectory
 
 meant for use with NMRlipids projects
 ------------------------------------------------------------
 Made by Joe,  Last edit 2017/02/02
------------------------------------------------------------
 input: Order parameter definitions
        gro and xtc file (or equivalents)
 output: order parameters (2 textfiles)
--------------------------------------------------------
"""

# coding: utf-8

import MDAnalysis as mda
import numpy as np
import math
import os, sys
import warnings
import re

from optparse import OptionParser
from collections import OrderedDict
from operator import add
from linecache import getline


# k_b = 0.0083144621  #kJ/Mol*K
# f_conc=55430  # factor for calculating concentrations of salts from numbers of ions/waters; in mM/L

bond_len_max = 1.5  # in A, max distance between atoms for reasonable OP calculation
bond_len_max_sq = bond_len_max**2


# %%
class OrderParameter:
    """
    :meta private:

    Class for storing&manipulating
    order parameter (OP) related metadata (definition, name, ...)
    and OP trajectories
    and methods to evaluate OPs.
    """

    def __init__(
        self, resname, atom_A_name, atom_B_name, M_atom_A_name, M_atom_B_name, *args
    ):  # removed name, resname from arguments
        """
        it doesn't matter which comes first,
        atom A or B, for OP calculation.
        """
        #        self.name = name             # name of the order parameter, a label
        self.resname = resname  # name of residue atoms are in
        self.atAname = atom_A_name
        self.atBname = atom_B_name
        self.M_atAname = M_atom_A_name
        self.M_atBname = M_atom_B_name
        self.name = (
            M_atom_A_name + " " + M_atom_B_name
        )  # generic name of atom A and atom B
        for field in self.__dict__:
            if not isinstance(field, str):
                warnings.warn(
                    "provided name >> {} << is not a string! \n \
                Unexpected behaviour might occur."
                )  # .format(field)
            else:
                if not field.strip():
                    raise RuntimeError(
                        "provided name >> {} << is empty! \n \
                    Cannot use empty names for atoms and OP definitions."
                    )  # .format(field)
        # extra optional arguments allow setting avg,std values -- suitable for reading-in results of this script
        if len(args) == 0:
            self.avg = None
            self.std = None
            self.stem = None
        elif len(args) == 2:
            self.avg = args[0]
            self.std = args[1]
            self.stem = None
        else:
            warnings.warn(
                "Number of optional positional arguments is {len}, not 2 or 0. Args: {args}\nWrong file format?"
            )
        self.traj = []  # for storing OPs
        self.selection = []

    def calc_OP(self, atoms):
        """
        calculates Order Parameter according to equation
        S = 1/2 * (3*cos(theta)^2 -1)
        """
        vec = atoms[1].position - atoms[0].position
        d2 = np.square(vec).sum()

        if d2 > bond_len_max_sq:
            at1 = atoms[0].name
            at2 = atoms[1].name
            resnr = atoms[0].resid
            d = math.sqrt(d2)
            warnings.warn(
                "Atomic distance for atoms \
            {at1} and {at2} in residue no. {resnr} is suspiciously \
            long: {d}!\nPBC removed???"
            )
        cos2 = vec[2] ** 2 / d2
        S = 0.5 * (3.0 * cos2 - 1.0)
        return S

    @property
    def get_avg_std_OP(self):
        """
        Provides average and stddev of all OPs in self.traj
        """
        # convert to numpy array
        return (np.mean(self.traj), np.std(self.traj))

    @property
    def get_avg_std_stem_OP(self):
        """
        Provides average and stddev of all OPs in self.traj
        """
        std = np.std(self.traj)
        # convert to numpy array
        return (np.mean(self.traj), std, std / np.sqrt(len(self.traj) - 1))

    @property
    def get_op_res(self):
        """
        Provides average and stddev of all OPs in self.traj
        """

        # convert to numpy array
        return self.traj


def read_trajs_calc_OPs(ordPars, top, trajs):
    """
    :meta private:

    procedure that
    creates MDAnalysis Universe with top,
    reads in trajectories trajs and then
    goes through every frame and
    evaluates each Order Parameter "S" from the list of OPs ordPars.
    ordPars : list of OrderParameter class
       each item in this list describes an Order parameter to be calculated in the trajectory
    top : str
        filename of a top file (e.g. conf.gro)
    trajs : list of strings
        filenames of trajectories
    """
    # read-in topology and trajectory
    mol = mda.Universe(top, trajs)

    # make atom selections for each OP and store it as its attribute for later use in trajectory
    for op in ordPars:
        # selection = pairs of atoms, split-by residues
        selection = mol.select_atoms(
            "resname {rnm} and name {atA} {atB}".format(
                rnm=op.resname, atA=op.atAname, atB=op.atBname
            )
        ).atoms.split("residue")
        #print(op.resname + " " + op.atAname + " " + op.atBname)

        for res in selection:
            # check if we have only 2 atoms (A & B) selected
            if res.n_atoms != 2:
                # print(res.resnames, res.resids)
                for atom in res.atoms:
                    # print(atom.name, atom.id)
                    atA = op.atAname
                    atB = op.atBname
                    nat = res.n_atoms
                    print(atA, atB, nat)
                    warnings.warn(
                        "Selection >> name {atA} {atB} << \
                   contains {nat} atoms, but should contain exactly 2!"
                    )
        op.selection = selection

    # go through trajectory frame-by-frame
    Nframes = len(mol.trajectory)
    for op in ordPars:
        #print(op.selection)
        Nres = len(op.selection)
        op.traj = [0] * Nres

    for frame in tqdm(mol.trajectory):
        for op in ordPars:  # .values():
            Nres = len(op.selection)
            for i in range(0, Nres):
                residue = op.selection[i]
                S = op.calc_OP(residue)
                op.traj[i] = op.traj[i] + S / Nframes


def parse_op_input(mapping_file, lipid_name):
    """:meta private:"""
    ordPars = []
    atomC = []
    atomH = []
    resname = lipid_name

    mapping_dict = {}
    with open("../BuildDatabank/mapping_files/" + mapping_file, "r") as yaml_file:
        mapping_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)
    yaml_file.close()

    regexp1_H = re.compile(r"M_[A-Z0-9]*C[0-9]*H[0-9]*_M")
    regexp2_H = re.compile(r"M_G[0-9]*H[0-9]*_M")
    regexp3_H = re.compile(r"M_C[0-9]*H[0-9]*_M")
    regexp1_C = re.compile(r"M_[A-Z0-9]*C[0-9]*_M")
    regexp2_C = re.compile(r"M_G[0-9]{1,2}_M")
    regexp3_C = re.compile(r"M_C[0-9]{1,2}_M")

    for mapping_key in mapping_dict.keys():
        #print(mapping_key)
        if (
            regexp1_C.search(mapping_key)
            or regexp2_C.search(mapping_key)
            or regexp3_C.search(mapping_key)
        ):
            atomC = [mapping_key, mapping_dict[mapping_key]["ATOMNAME"]]
            try:
                resname = mapping_dict[mapping_key]["RESIDUE"]
            except:
                pass
            atomH = []
        elif (
            regexp1_H.search(mapping_key)
            or regexp2_H.search(mapping_key)
            or regexp3_H.search(mapping_key)
        ):
            atomH = [mapping_key, mapping_dict[mapping_key]["ATOMNAME"]]
        else:
            atomC = []
            atomH = []

        #print(atomC)
        #print(atomH)

        if atomH and not len(atomC):
            print("Cannot define carbon for the hydrogen %s (%s)" % (atomH[0], atomH[1]), 
                file=sys.stderr )
            continue
        if atomH and len(atomC):
            #print(resname, atomC[1], atomH[1], atomC[0], atomH[0])
            items = [atomC[1], atomH[1], atomC[0], atomH[0]]
            op = OrderParameter(resname, items[0], items[1], items[2], items[3])
            #print(op.avg,op.resname,op.atAname,op.atBname,op.selection)
            ordPars.append(op)
    return ordPars


# def parse_op_input(fname):
#    """
#    parses input file with Order Parameter definitions
#    file format is as follows:
#    OP_name    resname    atom1    atom2
#    (flexible cols)
#    fname : string
#        input file name
#    returns : dictionary
#        with OrderParameters class instances
#    """
#    ordPars = OrderedDict()
#    try:
#        with open(fname,"r") as f:
#            for line in f.readlines():
#                if not line.startswith("#"):
#                    items = line.split()
#
#                    ordPars[items[0]] = OrderParameter(*items)
#
#    except:
#        inpf=opts.inp_fname
#        raise RuntimeError("Couldn't read input file >> {inpf} <<")
#    return ordPars


# %%


def find_OP(inp_fname, top_fname, traj_fname, lipid_name):
    """    :meta private:"""
    ordPars = parse_op_input(inp_fname, lipid_name)

    read_trajs_calc_OPs(ordPars, top_fname, traj_fname)

    return ordPars


#######################################################################################


def calc_angle(atoms, com):
    """
    :meta private:
    calculates the angle between the vector and z-axis in degrees
    no PBC check!
    Calculates the center of mass of the selected atoms to invert bottom leaflet vector
    """
    vec = atoms[1].position - atoms[0].position
    d = math.sqrt(np.square(vec).sum())
    cos = vec[2] / d
    # values for the bottom leaflet are inverted so that
    # they have the same nomenclature as the top leaflet
    cos *= math.copysign(1.0, atoms[0].position[2] - com)
    try:
        angle = math.degrees(math.acos(cos))
    except ValueError:
        if abs(cos) >= 1.0:
            print("Cosine is too large = {} --> truncating it to +/-1.0".format(cos))
            cos = math.copysign(1.0, cos)
            angle = math.degrees(math.acos(cos))
    return angle


def calc_z_dim(gro):
    """
    :meta private:
    Returns the simulation box dimension in z-direction from coordinate file.
    
    :param gro: coordinate in ``gro``, ``pdb`` or corresponding format.

    :return: size of box z-direction.
    """
    u = mda.Universe(gro)
    z = u.dimensions[2]
    return z

def system2MDanalysisUniverse(system):
    """
    Takes the ``system`` dictionary as an input, downloads the required files to the NMRlipids databank 
    directory and retuns MDAnalysis universe corressponding the ``system``.

    :param system: NMRlipids databank dictionary describing the simulation.

    :return: MDAnalysis universe
    """
    print(os.path.dirname(os.path.realpath(__file__)))
    systemPath = os.path.dirname(os.path.realpath(__file__)) + '/../../Data/Simulations/' + system['path'] 
    doi = system.get('DOI')
    skipDownloading: bool = (doi == 'localhost')
    if skipDownloading:
        print("NOTE: The system with 'localhost' DOI should be downloaded by the user.")

    trj = system.get('TRJ')
    trj_name = systemPath + system.get('TRJ')[0][0]
    software = system['SOFTWARE']

    if (skipDownloading):
        if (not os.path.isfile(trj_name)):
            raise FileNotFoundError(f"Trajectory should be downloaded [{trj_name}] by user")
    else:
        trj_url = resolve_download_file_url(doi, trj[0][0])                            
        if (not os.path.isfile(trj_name)):
            print('Downloading trajectory with the size of ', system['TRAJECTORY_SIZE'], ' to ', system['path'])
            response = urllib.request.urlretrieve(trj_url, trj_name)

    if 'gromacs' in software: 
        tpr = system.get('TPR')
        tpr_name = systemPath + tpr[0][0]
        if skipDownloading:
            if (not os.path.isfile(tpr_name)):
                raise FileNotFoundError(f"TPR should be downloaded [{tpr_name}]")
        else:
            tpr_url = resolve_download_file_url(doi, tpr[0][0])
            if (not os.path.isfile(tpr_name)):
                response = urllib.request.urlretrieve(tpr_url, tpr_name)
        try:
            u = mda.Universe(tpr_name, trj_name)
        except:
            conf = systemPath + '/conf.gro'
            if (not os.path.isfile(tpr_name)):
                print("Generating conf.gro because MDAnalysis cannot read tpr version")
                if 'WARNINGS' in system and 'GROMACS_VERSION' in system['WARNINGS'] and system['WARNINGS']['GROMACS_VERSION'] == 'gromacs3':
                    os.system('echo System | editconf -f '+ tpr_name + ' -o ' + conf)
                else:
                    os.system('echo System | gmx trjconv -s '+ tpr_name + ' -f '+ trj_name + ' -dump 0 -o ' + conf)
            u = mda.Universe(conf, trj_name)

    elif 'openMM' or 'NAMD' in software:
        pdb = system.get('PDB')
        pdb_name = systemPath + pdb[0][0]
        if skipDownloading:
            if (not os.path.isfile(pdb_name)):
                raise FileNotFoundError(f"PDB should be downloaded [{pdb_name}]")
        else:            
            pdb_url = resolve_download_file_url(doi, pdb[0][0])
            if (not os.path.isfile(pdb_name)):
                response = urllib.request.urlretrieve(pdb_url, pdb_name)
        u = mda.Universe(pdb_name, trj_name)

    else:
        raise NotImplementedError('Other than gromacs, openMM or NAMD are yet to be implemented.')
        
    return u


def read_trj_PN_angles(molname, atom1, atom2, MDAuniverse):
    """
    Calculates the P-N vector angles with respect to membrane normal from the simulation defined by the MDAnalysis universe. 

    :param molname: residue name of the molecule for which the P-N vector angle will be calculated
    :param atom1: name of the P atom in the simulation
    :param atom2: name of the N atom in the simulation
    :param MDAuniverse: MDAnalysis universe of the simulation to be analyzed

    :return: list where the first element are the angles of all molecules as a function of time, second element contains time averages for each molecule, third element contains the average angle over time and molecules, and fourth element is the error of the mean calculated over molecules. 
    """

    #systemPath = os.path.dirname(os.path.realpath(__file__)) + system['path'] 
    
    #if 'gromacs' in system['SOFTWARE']:
    #    print(systemPath, system['TPR'][0][0])
    #    top_fname = systemPath + system['TPR'][0][0]
    #else:
    #    top_fname = system['PDB']
    #traj_fname =  systemPath +  system['TRJ'][0][0]

    mol = MDAuniverse #mda.Universe(top_fname, traj_fname)

    selection = mol.select_atoms(
        "resname " + molname + " and (name " + atom1 + ")",
        "resname " + molname + " and (name " + atom2 + ")",
    ).atoms.split("residue")
    com = mol.select_atoms(
        "resname " + molname + " and (name " + atom1 + " or name " + atom2 + ")"
    ).center_of_mass()

    Nres = len(selection)
    Nframes = len(mol.trajectory)
    angles = np.zeros((Nres, Nframes))

    resAverageAngles = [0] * Nres
    resSTDerror = [0] * Nres
    j = 0

    for frame in mol.trajectory:
        for i in range(0, Nres):
            residue = selection[i]
            angles[i, j] = calc_angle(residue, com[2])
        j = j + 1
    for i in range(0, Nres):
        resAverageAngles[i] = sum(angles[i, :]) / Nframes
        resSTDerror[i] = np.std(angles[i, :])

    totalAverage = sum(resAverageAngles) / Nres
    totalSTDerror = np.std(resAverageAngles) / np.sqrt(Nres)

    return angles, resAverageAngles, totalAverage, totalSTDerror


###############################################################################################################


def calc_file_sha1_hash(fi: str, step: int = 4096) -> str:
    """
    :meta private:
    Calculates sha1 hash of given file using hashlib

    Args:
        fi (str): path to file
        step (int, optional): file read bytes step. Defaults to 4096.

    Returns:
        str: sha1 filehash of 40 char length
    """
    sha1_hash = hashlib.sha1()
    with open(fi, "rb") as f:
        with tqdm(total=math.ceil(os.path.getsize(fi) / step)) as pbar:
            # Read and update hash string value in blocks of 4K
            for byte_block in iter(lambda: f.read(step), b""):
                sha1_hash.update(byte_block)
                pbar.update(1)
    return sha1_hash.hexdigest()


def create_databank_directories(sim, sim_hashes, out) -> str:
    """
    :meta private:
    create nested output directory structure to save results

    Args:
        sim (_type_): Processed simulation entries
        sim_hashes (_type_): file hashes needed for directory structure
        out (str): output base path

    Raises:
        NotImplementedError: unsupported simulation software
        OSError: Error while creating the output directory

    Returns:
        str: output directory
    """
    # resolve output dir naming
    if sim["SOFTWARE"] == "gromacs":
        head_dir = sim_hashes.get("TPR")[0][1][0:3]
        sub_dir1 = sim_hashes.get("TPR")[0][1][3:6]
        sub_dir2 = sim_hashes.get("TPR")[0][1]
        sub_dir3 = sim_hashes.get("TRJ")[0][1]
    elif sim["SOFTWARE"] == "openMM" or sim["SOFTWARE"] == "NAMD":
        head_dir = sim_hashes.get("TRJ")[0][1][0:3]
        sub_dir1 = sim_hashes.get("TRJ")[0][1][3:6]
        sub_dir2 = sim_hashes.get("TRJ")[0][1]
        sub_dir3 = sim_hashes.get("TRJ")[0][1]
    else:
        raise NotImplementedError(f"sim software '{sim['SOFTWARE']}' not supported")

    directory_path = os.path.join(out, head_dir, sub_dir1, sub_dir2, sub_dir3)

    logger.debug(f"output_dir = {directory_path}")

    # destination directory is not empty
    if os.path.exists(directory_path) and os.listdir(directory_path) != 0:
        logger.warning(
            f"output directory '{directory_path}' is not empty. Data may be overriden."
        )

    # create directories
    os.makedirs(directory_path, exist_ok=True)

    return directory_path



class YamlBadConfigException(Exception):
    """
    :meta private:
    Custom Exception class for parsing the yaml configuration
    """

    def __init__(self, *args, **kwargs) -> None:
        Exception.__init__(self, *args, **kwargs)


def parse_valid_config_settings(info_yaml: dict) -> (dict, List[str]):
    """
    :meta private:
    Parses, validates and updates dict entries from yaml configuration file.

    Args:
        info_yaml (dict): info.yaml of database to add
    Raises:
        KeyError: Missing required key in info.yaml
        YamlBadConfigException: Incorrect or incompatible configuration
    Returns:
        dict: updated sim dict
        list[str]: list of filenames to download
    """

    sim = copy.deepcopy(info_yaml)  # mutable objects are called by reference in Python

    # STEP 1 - check supported simulation software
    if "SOFTWARE" not in sim:
        raise KeyError("'SOFTWARE' Parameter missing in yaml")

    if sim["SOFTWARE"].upper() in software_dict.keys():
        logger.info(f"Simulation uses supported software '{sim['SOFTWARE'].upper()}'")
    else:
        raise YamlBadConfigException(
            f"Simulation uses unsupported software '{sim['SOFTWARE'].upper()}'"
        )

    software_sim = software_dict[
        sim["SOFTWARE"].upper()
    ]  # related to dicts in this file

    # STEP 2 - check required keys defined by sim software used
    software_required_keys = [k for k, v in software_sim.items() if v["REQUIRED"]]

    # are ALL required keys are present in sim dict and defined (not of NoneType) ?
    if not all(
        (k in list(sim.keys())) and (sim[k] is not None) for k in software_required_keys
    ):
        missing_keys = [k for k in software_required_keys if k not in list(sim.keys()) ]
        raise YamlBadConfigException(
            f"Required '{sim['SOFTWARE'].upper()}' sim keys missing or not defined in conf file: {', '.join(missing_keys)}"
        )

    logger.debug(
        f"all {len(software_required_keys)} required '{sim['SOFTWARE'].upper()}' sim keys are present"
    )

    # STEP 3 - check working directory
    if "DIR_WRK" not in sim:
        raise KeyError("'DIR_WRK' Parameter missing in yaml")
    dir_wrk = sim["DIR_WRK"]

    # STEP 4 - Check that all entry keys provided for each simulation are valid
    files_tbd = []

    #   loop config entries
    for key_sim, value_sim in sim.items():
        logger.debug(f"processing entry: sim['{key_sim}'] = {str(value_sim)}")

        if key_sim.upper() in "SOFTWARE":  # skip 'SOFTWARE' entry
            continue

        # STEP 4.1.
        # Anne: check if key is in molecules_dict, molecule_numbers_dict or molecule_ff_dict too
        if (
            (key_sim.upper() not in software_sim.keys())
            and (key_sim.upper() not in molecules_dict.keys())
            and (key_sim.upper() not in lipids_dict.keys())
            and (key_sim.upper() not in molecule_ff_dict.keys())
        ):
            logger.error(
                f"key_sim '{key_sim}' in {sim['SOFTWARE'].lower()}_dict' : {key_sim.upper() in software_sim.keys()}"
            )
            logger.error(
                f"key_sim '{key_sim}' in molecules_dict : {key_sim.upper() in molecules_dict.keys()}"
            )
            logger.error(
                f"key_sim '{key_sim}' in lipids_dict : {key_sim.upper() in lipids_dict.keys()}"
            )
            logger.error(
                f"key_sim '{key_sim}' in molecule_ff_dict : {key_sim.upper() in molecule_ff_dict.keys()}"
            )
            raise YamlBadConfigException(
                f"'{key_sim}' not supported: Not found in '{sim['SOFTWARE'].lower()}_dict', 'molecules_dict', 'lipids_dict' and 'molecule_ff_dict'"
            )
        elif (
            key_sim.upper() not in software_sim.keys()
        ):  # hotfix for unkown yaml keys. TODO improve check 4.1?
            logger.warning(
                f"ignoring yaml entry '{key_sim}', not found in '{sim['SOFTWARE'].lower()}_dict'"
            )
            continue

        # STEP 4.2.
        # entries with files information to contain file names in arrays
        if "TYPE" in software_sim[key_sim]:
            if "file" in software_sim[key_sim]["TYPE"]:  # entry_type
                logger.debug(
                    f"-> found '{key_sim}:{software_sim[key_sim]}' of 'TYPE' file"
                )  # DEBUG

                if value_sim is None:
                    logger.debug(f"entry '{key_sim}' has NoneType value, skipping")
                # already a list -> ok
                elif isinstance(value_sim, list):
                    logger.debug(f"value_sim '{value_sim}' is already a list, skipping")
                    files_tbd.extend(value_sim)
                else:
                    value_sim_splitted = value_sim.split(";")

                    if len(value_sim_splitted) == 0:
                        raise YamlBadConfigException(
                            f"found no file to download for entry '{key_sim}:{software_sim[key_sim]}'"
                        )
                    # in case there are multiple files for one entry
                    elif len(value_sim_splitted) > 1:
                        files_list = []
                        for file_provided in value_sim.split(";"):
                            files_list.append([file_provided.strip()])
                        sim[
                            key_sim
                        ] = files_list  # replace ; separated string with list
                    else:
                        # print(f"value_sim_splitted = {value_sim_splitted}")
                        sim[key_sim] = [
                            [f.strip()] for f in value_sim_splitted
                        ]  # IMPORTANT: Needs to be list of lists for now
                    files_tbd.extend(f[0] for f in sim[key_sim])
                    # print(f"sim[{key_sim}] = {sim[key_sim]}")

                # STEP 4.3.
                # Batuhan: In conf file only one psf/tpr/pdb file allowed each (can coexist), multiple TRJ files are ok
                # TODO true for all sim software?
                # TODO add dict entry param "unique" instead?
                if key_sim.upper() in ["PSF", "TPR", "PDB"] and len(sim[key_sim]) > 1:
                    raise YamlBadConfigException(
                        f"only one '{key_sim}' entry file allowed, but got {len(sim[key_sim])}: {sim[key_sim]}"
                    )

        else:
            logger.warn(
                f"skipping key '{key_sim}': Not defined in software_sim library"
            )

    logger.info(
        f"found {len(files_tbd)} ressources to download: {', '.join(files_tbd)}"
    )

    return sim, files_tbd

def calcArea(system):
    """
    Returns area of the calculated based on the area per lipid stored in the databank.

    :param system: a system dictionary

    :return: area of the system (Å^2)
    """
    APL = CalcAreaPerMolecule(system)
    Nlipid = 0
    for molecule in system['COMPOSITION']:
        if molecule in lipids_dict:
            Nlipid += np.sum(system['COMPOSITION'][molecule]['COUNT'])
    print(Nlipid, APL)
    return Nlipid*APL/2

def GetFormFactorMin(system):
    """
    Return list of minima of form factor of ``system``.

    :param system: a system dictionary

    :return: list of form factor minima
    """
    FormFactorPath = os.path.dirname(os.path.realpath(__file__)) + '/../../Data/Simulations/' + system['path'] + 'FormFactor.json'
    #try:
    f = open(FormFactorPath)
    FormFactor = json.load(f)
    min = 1000
    iprev = FormFactor[0][1]
    iprevD = 0
    minX = []
    for i in FormFactor:
        iD = i[1]-iprev
        if iD > 0 and iprevD < 0 and i[0] > 0.1:
            minX.append(i[0])
        iprevD = i[1]-iprev
        iprev = i[1]
        
    return(minX)


def averageOrderParameters(system):
    """
    Returns average order paramaters of sn-1 and sn-2 acyl chains based on universal atom names. The names starting with M_G1C will be assigned to sn-1 and names starting M_G2C to sn-2. 

    :parameters system: a system dictionary

    :return: average of sn-1 and sn-2 order parameters
    """
    
    #DataBankPath = '../../Databank/Data/'
    path = os.path.dirname(os.path.realpath(__file__)) + '/../../Data/Simulations/' + system['path']
    #DataBankPath + 'Simulations/' +system['path']
    
    sn1sum = 0
    sn1count = 0
    sn2sum = 0
    sn2count = 0
    
    for lipid in system['COMPOSITION']:
        if lipid in lipids_dict and not 'CHOL' in lipid:
            OPpathSIM = path + lipid + 'OrderParameters.json'
            with open(OPpathSIM) as json_file:
                OPsim = json.load(json_file)
    
            for key in OPsim:
                if 'M_G1C' in key:
                    sn1sum += float(OPsim[key][0][0])
                    sn1count += 1
                elif 'M_G2C' in key:
                    sn2sum += float(OPsim[key][0][0])
                    sn2count += 1
                    
    return sn1sum/sn1count, sn2sum/sn2count

def calcLipidFraction(system, lipid):
    """
    Returns the number fraction of ``lipid`` with respect to total number of lipids.

    :param system: a system dictionary
    :param lipid: universal molecule name of lipid

    :return: number fraction of ``lipid`` with respect total number of lipids
    """
    NlipidTOT = 0
    for molecule in system['COMPOSITION']:
        if molecule in lipids_dict:
            NlipidTOT += np.sum(system['COMPOSITION'][molecule]['COUNT'])
    
    Nlipid = 0
    for molecule in system['COMPOSITION']:
        if lipid in molecule:
            Nlipid += np.sum(system['COMPOSITION'][molecule]['COUNT'])
            
    return Nlipid/NlipidTOT


def getHydrationLevel(system):
    """
    Returns hydration level of the system, i.e., number of water molecules divided by number of lipid molecules.

    :param system: a system dictionary

    :return: number of water molecules divided by number of lipid molecules
    """
    Nlipid = 0
    for molecule in system['COMPOSITION']:
        if molecule in lipids_dict:
            Nlipid += np.sum(system['COMPOSITION'][molecule]['COUNT'])
    Nwater = system['COMPOSITION']['SOL']['COUNT']
    return Nwater/Nlipid
