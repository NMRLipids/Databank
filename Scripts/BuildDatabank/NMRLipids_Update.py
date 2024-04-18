#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 14:27:05 2021

@author: Fabs

"""

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODULES
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

import os
import re
import sys
import glob
import json
import yaml
import pymysql
import argparse
import numpy as np


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ARGUMENTS
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Program description
parser = argparse.ArgumentParser(description =
    '''NMRLipids Update v0.1''')

# Ubication of data
parser.add_argument( "-f", "--folder", type= str, default = 'Databank-main',
    help= ''' Absolute path of the copy of the repository.\n
    Default: %(default)s ''' )

parser.add_argument( "-sf", "--simulation_folder", type = str, default = 'Data/Simulations/',
    help = ''' Path of the simulation folder in the repositoty.\n
    Default: %(default)s ''' )
parser.add_argument( "-rf", "--ranking_folder", type = str, default = 'Data/Ranking/',
    help = ''' Path of the ranking folder in the repositoty.\n
    Default: %(default)s ''' )
parser.add_argument( "-ef", "--experiment_folder", type = str, default = 'Data/experiments/',
    help = ''' Path of the experiments folder in the repositoty.\n
    Default: %(default)s ''' )

parser.add_argument( "-c", "--config", type=str, default = "config.json",
    help = ''' JSON file with the configuration of the connection to the DB.
    Default: %(default)s ''')

# System properties
parser.add_argument( "-s", "--systems", type = str, nargs = '+', # REQUIRED
    help = """ Path of the system(s). """ )

args = parser.parse_args()


sys.path.insert(1, args.folder + '/Scripts/BuildDatabank/')
import databank_defs as NMRDict

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#  Elements from databankLibrary
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SQL Queries
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def SQL_Select( Table: str, Values: list, Condition: dict = {} ) -> str:
    '''
    Generate a SQL query to select values in a table

    Parameters
    ----------
    Table : str
        Name of the table.
    Values : list
        List of values to select.
    Condition : dict, optional
        Condition(s) for the search. 

    Returns
    -------
    str
        The SQL query:
        SELECT Values[0], (...), Values[-1] FROM Table 
          WHEN Condition.keys()[0]=Condition.value()[0] AND ... Condition.keys()[-1]=Condition.value()[-1]
    '''
    String = ' SELECT `' + '`, `'.join( Values ) + '` FROM `' + Table + '`'
    
    # Add a condition to the search
    if Condition:  String += ' WHERE `' + '" AND `'.join( [ str( Par ) + '` = "' + str( Condition[Par] ) for Par in Condition ] ) + '"'

    return String


def SQL_Create( Table: str, Values: dict, Condition: dict = {} ) -> str:
    '''
    Generate a SQL query to insert a new entry in a table.

    Parameters
    ----------
    Table : str
        Name of the table.
    Values : dict
        List of values to insert.
    Condition : dict, optional
        Condition(s) for the insertion. 

    Returns
    -------
    str
        The SQL query:
        INSERT INTO Table ( Values.keys()[0], ..., Values.keys()[-1] ) VALUES ( Values.values()[0], ..., Values.values()[-1] ) 
          WHEN Condition.keys()[0]=Condition.value()[0] AND ... Condition.keys()[-1]=Condition.value()[-1]
    '''
    String = ' INSERT INTO `' + Table + '` (`' +'`, `'.join( [ Par for Par in Values ] ) + '`) VALUES ("' + '", "'.join( [ str( Values[Par] ) for Par in Values ] ) + '")'
    
    # Add a condition to the search
    if Condition:  String += ' WHERE `' + '" AND `'.join( [ str( Par ) + '` = "' + str( Condition[Par] ) for Par in Condition ] ) + '"'

    return String


def SQL_Update( Table: str, Values: dict, Condition: dict = {} ) -> str:
    '''
    Generate a SQL query to update an entry in a table.

    Parameters
    ----------
    Table : str
        Name of the table.
    Values : dict
        List of values to insert.
    Condition : dict, optional
        Condition(s) for the insertion. 

    Returns
    -------
    str
        The SQL query:
        UPDATE Table SET Values.keys()[0] = Values.values()[0], ..., Values.keys()[-1] = Values.values()[-1]
          WHEN Condition.keys()[0]=Condition.value()[0] AND ... Condition.keys()[-1]=Condition.value()[-1]
    '''
    String = ' UPDATE `' + Table + '` SET `' + '", `'.join( [ str( Par ) + '` = "' + str( Values[Par] ) for Par in Values ] ) + '"'
    
    # Add a condition to the search
    if Condition:  String += ' WHERE `' + '" AND `'.join( [ str( Par ) + '` = "' + str( Condition[Par] ) for Par in Condition ] ) + '"'

    return String


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Functions
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def CheckEntry( Table: str, Information: dict = {} ) -> int:
    '''
    Find an entry in the DB

    Parameters
    ----------
    Table : str
        Name of the table.
    Information : dict, optional
        Values to check.

    Returns
    -------
    int or None
        ID of the entry in the table. If it does not exists, the value is None
    '''
    
    # Create a cursor
    cursor = database.cursor()
    
    # Find the ID(s) of the entry matching the condition
    cursor.execute( SQL_Select( Table, [ "id" ], Information ) )
    
    # The list of IDs
    ID = cursor.fetchall()
    
    # More than ID will raise an error
    if len(ID)!=1:
        return None
    # If only 1 ID, extract its value
    else:
        if type(ID[0]) == int:
            return ID[0]
        else: 
            return ID[0][0]


def CreateEntry( Table: str, Information: dict ) -> int:
    '''
    Add an entry into a table

    Parameters
    ----------
    Table : str
        Name of the table.
    Information : dict, optional
        Values to add.

    Returns
    -------
    int
        ID of the entry in the table. If it does not work, value will be 0.
    '''
    
    # Function
    #   CreateEntry
    # Add a new entry to a table in the DB.
    
    # Create a cursor
    cursor = database.cursor()
    
    # Execute the query creating a new entry
    cursor.execute( SQL_Create( Table, Information ) )
    
    # Commit the changes
    database.commit()
    
    # Check if the entry was created
    ID = CheckEntry( Table, Information )
    
    # If there is not an ID, raise an error (the table was not created)
    if not ID:
        print("WARNING: Something may have gone wrong with the table {}".format(Table))
        return 0
    # If an ID is obtained, the entry was created succesfuly
    else:
        print("A new entry was created in {}: index {}".format(Table, ID))
        return ID


def UpdateEntry( Table: str, Information: dict, Condition: dict ):
    '''
    Updates an entry in a table.

    Parameters
    ----------
    Table : str
        Name of the table.
    Information : dict
        Values to add.
    Condition : dict
        Conditions to select the entry.
    '''
    
    # Create a cursor
    cursor = database.cursor()
    
    # Execute the query updating an entry
    cursor.execute( SQL_Update( Table, Information, Condition ) )
    
    # Commit the changes
    database.commit()
    
    return print("Entry {} in table {} was updated".format(Condition["id"], Table))
    
    
def DBEntry( Table: str, Information: dict, Minimal: dict = {} ) -> tuple:
    '''
    Manages entries in the DB. If the Minimal information is not found in an
    existing entry, a new one is created; else, the matching entry is updated
    when a discrepancy between the minimal and the total information appears.

    Parameters
    ----------
    Table : str
        Name of the table
    Information : dict
        Total infomation of the entry.
    Minimal : dict, optional
        Minimal information of the entry.

    Returns
    -------
    tuple
        The ID of the created/updated entry.
    '''
    
##### TEMPORAL ##### 
    # Delete the Nan from the Information dictionary
    # The presence of NaN in the DB leads to problems dealing with the data
    # A solution must be found for this problem in the future.
    for entry in Information:
        if Information[ entry ] == "nan" : Information[ entry ] = 0
####################
    
    # Check the existence of the entry
    EntryID = CheckEntry( Table, Minimal )
    
    # If there is an ID associated to the minimal information...
    if EntryID:
        
        # Check if the Minimal information matches the whole Information
        # If they don't match...
        if Information != Minimal:
            
            # Add the ID of the entry to the minimal information
            Minimal["id"] = EntryID
            
            # Update the entry
            UpdateEntry( Table, Information, Minimal )
            return EntryID
        
        # If the minimal information and the total information match, no update is required.
        else:
            print("It was not necessary to update entry {} in table {}".format( EntryID, Table ) )
            return EntryID
    
    # If there is not an entry...
    else:
        # Create a new one with the data from Information
        EntryID = CreateEntry( Table, Information )
        return EntryID
    

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MAIN PROGRAM
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

HETEROMOLECULES_LIST = [ "CHOL", "DCHOL", "C20", "C30" ]
FAILS = []

if __name__ == '__main__':
    # Load the configuration of the connection
    config = json.load( open( args.config, "r" ) )
    
    # Conect to the database
    database = pymysql.connect( **config )
    
    # Format of the folder
    args.folder = os.path.abspath( args.folder ) + '/'
    
    # Path to the ranking files
    PATH_RANKING = args.folder + args.ranking_folder
    
    # If -s(ystems) is given, only modify that entries
    try:
        systems = []
        for system in args.systems:
            # Format the system info
            if system[0] == '/':
                system = system[1:]
            if system[-1] != '/':
                system += '/'
            
            # The path to the simulation file
            PATH_SIMULATION = args.folder + args.simulation_folder + system
                
            # Read the README.yaml file, the information of the simulation
            with open( PATH_SIMULATION + 'README.yaml' ) as File:
                README = yaml.load( File, Loader = yaml.FullLoader )
                
            # Add the math of the file
            README["path"] = system    
                
            systems.append( README )
        
    # Otherwise load the whole databank
    except:
        systems = initialize_databank( args.folder )    

    # Iterate over the loaded systems
    for README in systems:
        
        try:
            print("\nCollecting data from system:")
            print( README["path"] + "\n" )
            
            # The location of the files
            PATH_SIMULATION = args.folder + args.simulation_folder + README["path"]    
        
            # In the case a field in the README does not exist, set its value to 0
            for field in ['AUTHORS_CONTACT','COMPOSITION','CPT','DATEOFRUNNING','DIR_WRK',
                          'DOI','EXPERIMENT','FF','FF_DATE','FF_SOURCE','GRO','LOG','NUMBER_OF_ATOMS',
                          'PREEQTIME','PUBLICATION','SOFTWARE','SOFTWARE_VERSION','SYSTEM','TEMPERATURE',
                          'TIMELEFTOUT','TOP','TPR','TRAJECTORY_SIZE','TRJ','TRJLENGTH','TYPEOFSYSTEM','WARNINGS','ID']:
                if not field in README: README[field] = 0        
            
### TABLE 1: forcefields
            # Collect the information about the forcefield
            Info = { "name":   README["FF"],
                     "date":   README["FF_DATE"], 
                     "source": README["FF_SOURCE"] }
            
            # Entry in the DB with the info of the FF
            FF_ID = DBEntry( 'forcefields', Info, Info )
        
        
### TABLE 2: lipids
            # Empty dictionaries for the info of the lipids
            Lipids = {}; Lipids_ID = {}; Lipid_Ranking = {}; Lipid_Quality = {}
            # Find the lipids in the composition
            for key in README["COMPOSITION"]:        
                if key in NMRDict.lipids_dict and not key in HETEROMOLECULES_LIST:
                    # Save the quality of the lipid
                    Store = True
                    
                    # Collect the info of the lipids
                    Info = { "forcefield_id": FF_ID,
                             "molecule":      key,
                             "name":          README["COMPOSITION"][ key ]["NAME"],
                             "mapping":       README["COMPOSITION"][ key ]["MAPPING"] }
        
                    # Entry in the DB with the info of the lipid
                    Lip_ID = DBEntry( 'lipids', Info, Info )
                    
                    # Store information for further steps
                    Lipids[ key ] = README["COMPOSITION"][ key ]["COUNT"]
                    Lipids_ID[ key ] = Lip_ID
        
##### TEMPORAL ##### 
# Must be chenged when the final structure is ready
# If the LIPID_FragmentQuality.json file will be defined for every system
# this part can be deleted and just read the quality (try at the end of this
# part). At this moment the ranking is not necessary in the DB, the web 
# already provides the result sorted by quality.
# Same for the heteromolecules.
                    Lipid_Ranking[ key ] = {}
                    # Find the position of the system in the ranking
                    for file in glob.glob( PATH_RANKING + key + "*" ):
        
                        #print(file)
                        
                        # Kind of ranking (total, headgroup...)
                        kind = re.search( '_(.*)_', file ).group(1)
                        
                        # Open the ranking file
                        with open( file ) as FILE:
                            RANKING_LIST = json.load( FILE )
                            
                            # Find the position of the system in the ranking
                            for SIM in range(len(RANKING_LIST)):
        
                                if README["path"] in RANKING_LIST[SIM]["system"]["path"]:
                                    Lipid_Ranking[ key ][kind] = SIM + 1
                                    
                                    if Store:
                                        Lipid_Quality[ key ] = RANKING_LIST[SIM][key]
                                        Store=False
                            
                            # If it does not have an assigned value, use None
                            if not kind in Lipid_Ranking[ key ]:
                                Lipid_Ranking[ key ][kind] =  0 
                        
                        FILE.close()
                      
                    # Read the quality file of the lipid (this will remain if the rest is removed)
                    try:
                        with open( PATH_SIMULATION + key + '_FragmentQuality.json' ) as FILE:
                            Lipid_Quality[ key ] = json.load( FILE )
                        FILE.close()
                    except:
                        Lipid_Quality[key] = { "total": 0,
                                               "headgroup": 0,
                                               "sn-1": 0,
                                               "sn-2": 0 }
                        
                    for t in ["total", "headgroup", "sn-1", "sn-2"]:
                        try:    Lipid_Quality[key][t] = Lipid_Quality[key][t]  if not np.isnan( Lipid_Quality[key][t] ) else 0
                        except: Lipid_Quality[key][t] = 0
####################
                        
    
### TABLE 3: ions
            # Empty dictionary for the info of the ions
            Ions = {}
            # Find the ions in the composition
            for key in README["COMPOSITION"]:
                if key in NMRDict.molecules_dict and key != "SOL" and not key in HETEROMOLECULES_LIST:
                    
                    # Collect the info of the ions
                    Info = { "forcefield_id": FF_ID,
                             "molecule":      key,
                             "name":          README["COMPOSITION"][ key ]["NAME"],
                             "mapping":       README["COMPOSITION"][ key ]["MAPPING"] }
                    
                    # Entry in the DB with the info of the ion
                    Ion_ID = DBEntry( 'ions', Info, Info )
                    
                    # Store information for further steps: Ions[name]=[ID,number]
                    Ions[ key ] = [ Ion_ID, README["COMPOSITION"][ key ]["COUNT"] ]
        
            
        
### TABLE 4: water_models (rename to solvent better?)
            if "SOL" in README["COMPOSITION"]:
                WatName = README["COMPOSITION"]["SOL"]["NAME"]
                
                WatNum = README["COMPOSITION"]["SOL"]["COUNT"]
            
                # Collect the info on the water
                Info = { "short_name":    WatName,
                         "mapping":       README["COMPOSITION"]["SOL"]["MAPPING"] }
                
                # Get an entry in the DB with the info of the water
                Wat_ID = DBEntry( 'water_models', Info, Info )
            
        
### TABLE 5: heteromolecules
# Heteromolecules are defined as the lipids for whom a distinction between the
# different parts (headgroup, sn-1, sn-2) was not made. They could be included
# in the lipids table and leaving some fields empty.
    
            # Empty dictionaries for the info of the heteromolecules
            Heteromolecules = {}; Heteromolecules_ID = {}; Heteromolecules_Ranking = {} 
            
            # Form factor
            Store = True; Heteromolecules_Quality = {}
        
            # Find the heteromolecules in the composition
            for key in README["COMPOSITION"]:        
                if key in NMRDict.lipids_dict and key in HETEROMOLECULES_LIST:
                    
                    # Collect the info of the lipids
                    Info = { "forcefield_id": FF_ID,
                             "molecule":      key,
                             "name":          README["COMPOSITION"][ key ]["NAME"],
                             "mapping":       README["COMPOSITION"][ key ]["MAPPING"] }
                    
                    # Entry in the DB with the info of the heteromolecules
                    Mol_ID = DBEntry( 'heteromolecules', Info, Info )
                    
                    # Store information for further steps
                    Heteromolecules[ key ] = README["COMPOSITION"][ key ]["COUNT"]
                    Heteromolecules_ID[ key ] = Mol_ID
        
##### TEMPORAL ##### 
# See the lipids table for the reasons
                    Heteromolecules_Ranking[ key ] = {}
                    
                    # Find the position of the system in the raking
                    for file in glob.glob( PATH_RANKING + key + "*" ):
                        
                        # Type of ranking
                        kind = re.search( '_(.*)_', file ).group(1)
                        
                        # Open the ranking file
                        with open( file ) as FILE:
                            RANKING_LIST = json.load( FILE )
                            
                            # Find the position of the system in the ranking
                            for SIM in range(len(RANKING_LIST)):
                                if README["path"] in RANKING_LIST[SIM]["system"]["path"]:
                                    Heteromolecules_Ranking[ key ][kind] = SIM + 1
                                    
                                    if Store:
                                        Heteromolecules_Quality[ key ] = RANKING_LIST[SIM][key]
                                        Store=False
                            
                            # If it does not have an assigned value, use None
                            if not kind in Heteromolecules_Ranking[ key ]:
                                Heteromolecules_Ranking[ key ][kind] =  0
                        
                        FILE.close()
                        
                    try:
                        with open( PATH_SIMULATION + key + '_FragmentQuality.json' ) as FILE:
                            Heteromolecules_Quality[ key ] = json.load( FILE )
                        FILE.close()
                    except:
                        Heteromolecules_Quality[key] = { "total": 0,
                                                         "headgroup": 0,
                                                         "tail": 0 }
                        
                    for t in ["total","headgroup","tail"]:
                        try:    Heteromolecules_Quality[key][t] = Heteromolecules_Quality[key][t] if not np.isnan( Heteromolecules_Quality[key][t] ) else 0
                        except: Heteromolecules_Quality[key][t] = 0
####################
        
        
### TABLE 6: membranes
            # Find the proportion of each lipid in the leaflets
            Names = [ [], [] ]; Number = [ [], [] ]
            
            for lipid in Lipids:
                if Lipids[lipid][0]: Names[0].append( lipid ); Number[0].append( str( Lipids[lipid][0] ) )
                if Lipids[lipid][1]: Names[1].append( lipid ); Number[1].append( str( Lipids[lipid][1] ) )
                
            for hetero in Heteromolecules:
                if Heteromolecules[hetero][0]: Names[0].append( hetero ); Number[0].append( str( Heteromolecules[hetero][0] ) )
                if Heteromolecules[hetero][1]: Names[1].append( hetero ); Number[1].append( str( Heteromolecules[hetero][1] ) ) 
            
            Names = [ ':'.join( Names[0] ), ':'.join( Names[1] ) ]
            Number = [ ':'.join( Number[0] ), ':'.join( Number[1] ) ]
            
            # Collect the information about the membrane
            Info = { "forcefield_id":   FF_ID, 
                     "lipid_names_l1":  Names[0],
                     "lipid_names_l2":  Names[1], 
                     "lipid_number_l1": Number[0], 
                     "lipid_number_l2": Number[1],
                     "geometry":        README["TYPEOFSYSTEM"] }
            
            # Entry in the DB with the info of the membrane
            Mem_ID = DBEntry( 'membranes', Info, Info )
            
            
### TABLE 7: trajectories
            # Collect the information about the simulation
            Info = { "id":              README["ID"],
                     "forcefield_id":   FF_ID, 
                     "membrane_id":     Mem_ID, 
                     "git_path":        README["path"],
                     "system":          README["SYSTEM"],
                     "author":          README["AUTHORS_CONTACT"],
                     "date":            README["DATEOFRUNNING"],
                     "dir_wrk":         README["DIR_WRK"],
                     "doi":             README["DOI"],
                     "number_of_atoms": README["NUMBER_OF_ATOMS"],
                     "preeq_time":      README["PREEQTIME"],
                     "publication":     README["PUBLICATION"],
                     "software":        README["SOFTWARE"],
                     "temperature":     README["TEMPERATURE"],
                     "timeleftout":     README["TIMELEFTOUT"],
                     "trj_size":        README["TRAJECTORY_SIZE"],
                     "trj_length":      README["TRJLENGTH"] }
        
            # The information that defines the trajectory
            Minimal = { "id":            README["ID"],
                        "forcefield_id": FF_ID, 
                        "membrane_id":   Mem_ID, 
                        "git_path":      README["path"],
                        "system":        README["SYSTEM"] }
            
            # Entry in the DB with the info of the trajectory
            Trj_ID = DBEntry( 'trajectories', Info, Minimal )    
        
            
### TABLE 8: trajectories_lipids
            TrjL_ID = {}
            for lipid in Lipids:
                # Collect the information of each lipid in the simulation
                Info = { "trajectory_id": Trj_ID, 
                         "lipid_id":      Lipids_ID[ lipid ],
                         "lipid_name":    lipid, 
                         "leaflet_1":     Lipids[ lipid ][0], 
                         "leaflet_2":     Lipids[ lipid ][1] }
                
                # The minimal information that identifies the lipid
                Minimal = { "trajectory_id": Trj_ID, 
                            "lipid_id":      Lipids_ID[ lipid ] }
                  
                
                # Entry in the DB with the info of the lipids in the simulation
                TrjL_ID[ lipid ] = DBEntry( 'trajectories_lipids', Info, Minimal )
            
            
### TABLE 9: trajectories_water
            if "SOL" in README["COMPOSITION"]:
                # Collect the information of the water in the simulation
                Info = { "trajectory_id": Trj_ID, 
                         "water_id":      Wat_ID, 
                         "water_name":    WatName, 
                         "number":        WatNum }
                
                # The minimal information of the water in the simulation
                Minimal = { "trajectory_id": Trj_ID, 
                            "water_id":      Wat_ID }
                
                # Entry in the DB with the info of the water in the simulation
                _ = DBEntry( 'trajectories_water', Info, Minimal )
            
        
### TABLE 10: trajectories_ions
            TrjI_ID = {}
            for ion in Ions:
                # Collect the information of each ion in the simulation
                Info = { "trajectory_id": Trj_ID, 
                         "ion_id":        Ions[ ion ][0],
                         "ion_name":      ion, 
                         "number":        Ions[ ion ][1] }
                
                # The minimal information that identifies the ion
                Minimal = { "trajectory_id": Trj_ID, 
                            "ion_id":        Ions[ ion ][0] }
                
                # Entry in the DB with the info of the ions in the simulation
                TrjI_ID[ ion ] = DBEntry( 'trajectories_ions', Info, Minimal )
        
        
### TABLE 11: trajectories_heteromolecules
            TrjM_ID = {}
            for hetero in Heteromolecules:
                # Collect the information of each heteromolecule in the simulation
                Info = { "trajectory_id": Trj_ID, 
                         "molecule_id":   Heteromolecules_ID[ hetero ],
                         "molecule_name": hetero, 
                         "leaflet_1":     Heteromolecules[ hetero ][0], 
                         "leaflet_2":     Heteromolecules[ hetero ][1] }
                
                # The minimal information that identifies the heteromolecule
                Minimal = { "trajectory_id": Trj_ID, 
                            "molecule_id":   Heteromolecules_ID[ hetero ] }
                
                # Entry in the DB with the info of the heteromolecules in the simulation
                TrjM_ID[ hetero ] = DBEntry( 'trajectories_heteromolecules', Info, Minimal )
            
        
### TABLE 12: trajectories_analysis
            # Find the bilayer thickness
            try:
                with open( PATH_SIMULATION + 'thickness.json' ) as FILE:
                    BLT = json.load( FILE )
                FILE.close()
            except:
                BLT = 0
            
            # Find the area per lipid
            try:
                with open( PATH_SIMULATION + 'apl.json' ) as FILE:
                    # Load the file
                    ApL = json.load( FILE )
                    
                    # Transform the dictionary into an array
                    ApL = np.array( [ [ float(key), float( ApL[key] ) ] for key in ApL ] )
                    
                    # Perform the mean
                    APL = np.mean( ApL[int(len(ApL[:,0])/2):,1] )
            except:
                APL = 0
            
            # Form factor quality
            try:
                with open( PATH_SIMULATION + 'FormFactorQuality.json' ) as FILE:
                    FFQ = json.load( FILE )
                FILE.close()
            except:
                FFQ = [4242, 0]
            
            # Read the quality file for the whole system
            try: 
                with open( PATH_SIMULATION + 'SYSTEM_quality.json' ) as FILE:
                    QUALITY_SYSTEM = json.load( FILE )
                FILE.close()
            except:
                QUALITY_SYSTEM = { "total": 0,
                                   "headgroup": 0,
                                   "tails": 0 }
            
            try:    FFExp = args.experiment_folder + "FormFactors/" + README["EXPERIMENT"]["FORMFACTOR"]
            except: FFExp = ''
            
            # Collect the information of the analysis of the trajectory
            Info = { "trajectory_id":          Trj_ID,
                     "bilayer_thickness":      BLT,
                     "area_per_lipid":         APL,
                     "area_per_lipid_file":    args.simulation_folder + README["path"] + 'apl.json',
                     "form_factor_file":       args.simulation_folder + README["path"] + 'FormFactor.json',
                     "quality_total":          QUALITY_SYSTEM["total"],
                     "quality_headgroups":     QUALITY_SYSTEM["headgroup"],
                     "quality_tails":          QUALITY_SYSTEM["tails"],
                     "form_factor_experiment": FFExp,
                     "form_factor_quality":    FFQ[0],
                     "form_factor_scaling":    FFQ[1] }
            
            # Collect the minimal information of the analysis of the trajectory
            Minimal = { "trajectory_id": Trj_ID } 
            
            # Entry in the DB with the info of the analysis of the simulation
            _ = DBEntry( 'trajectories_analysis', Info, Minimal )
            
        
### TABLE 13: trajectories_analysis_lipids
            for lipid in Lipids:
                try:    OPExp = args.experiment_folder + "OrderParameters/" + [ i for i in README["EXPERIMENT"]["ORDERPARAMETER"][lipid].values()][0] + '/' + lipid + '_Order_Parameters.json'
                except: OPExp = ''
                
                # Collect the information of each lipid in the simulation
                Info = { "trajectory_id":                Trj_ID,
                         "lipid_id":                     Lipids_ID[ lipid ],
                         "quality_total":                Lipid_Quality[ lipid ]["total"] ,
                         "quality_hg":                   Lipid_Quality[ lipid ]["headgroup"],
                         "quality_sn-1":                 Lipid_Quality[ lipid ]["sn-1"],
                         "quality_sn-2":                 Lipid_Quality[ lipid ]["sn-1"],
                         "order_parameters_file":        args.simulation_folder + README["path"] + lipid + 'OrderParameters.json',
                         "order_parameters_experiment":  OPExp,
                         "order_parameters_quality":     args.simulation_folder + README["path"] + lipid + '_OrderParameters_quality.json' }
                
                # The minimal information that identifies the lipid in the simulation
                Minimal = { "trajectory_id": Trj_ID,
                            "lipid_id":      Lipids_ID[ lipid ] }
                
                # Entry in the DB with the info of the analysis of the lipid in the simulation
                _ = DBEntry( 'trajectories_analysis_lipids', Info, Minimal )
            
        
### TABLE 14: trajectories_analysis_heteromolecules
            for hetero in Heteromolecules:
                try:    OPExp = args.experiment_folder + "OrderParameters/" + [ i for i in README["EXPERIMENT"]["ORDERPARAMETER"][hetero].values()][0] + '/' + hetero + '_Order_Parameters.json'
                except: OPExp = ''       
        
                # Collect the information of each heteromolecule in the simulation
                Info = { "trajectory_id":                Trj_ID,
                         "molecule_id":                  Heteromolecules_ID[ hetero ],
                         "quality_total":                Heteromolecules_Quality[ hetero ]["total"],
                         "quality_hg":                   Heteromolecules_Quality[ hetero ]["headgroup"],
                         "quality_tails":                Heteromolecules_Quality[ hetero ]["tail"],
                         "order_parameters_file":        args.simulation_folder + README["path"] + hetero + 'OrderParameters.json' ,
                         "order_parameters_experiment":  OPExp,
                         "order_parameters_quality":     args.simulation_folder + README["path"] + hetero + '_OrderParameters_quality.json' }
        
                # The minimal information that identifies the heteromolecule in the simulation
                Minimal = { "trajectory_id": Trj_ID,
                            "molecule_id":   Heteromolecules_ID[ hetero ] }
                
                # Entry in the DB with the info of the analysis of the heteromolecule in the simulation
                _ = DBEntry( 'trajectories_analysis_heteromolecules', Info, Minimal )
        
        
### TABLE 15: trajectory_analysis_ions
            for ion in Ions:
                # Collect the information of the ions in the simulation
                Info = { "trajectory_id": Trj_ID,
                         "ion_id":        Ions[ ion ][0] }
        
                # The minimal information that identifies the ion in the simulation
                #Minimal = { "trajectory_id": Trj_ID,
                #            "ion_id":        Ions[ ion ][0] }
                
                # Entry in the DB with the info of the analysis of the ion in the simulation
                _ = DBEntry( 'trajectories_analysis_ions', Info, Info )
                
        
### TABLE 16: trajectory_analysis_water
            if "SOL" in README["COMPOSITION"]:
                # Collect the information of the water in the simulation
                Info = { "trajectory_id": Trj_ID,
                         "water_id":      Wat_ID, }
                         # "density_file":  args.densities_folder + README["path"] + '/SOLdensity.xvg' }
                
                # Entry in the DB with the info of the analysis of the water in the simulation
                _ = DBEntry( 'trajectories_analysis_water', Info, Info )
        
        
##### TEMPORAL ##### 
# The table may be removed in the future, and the quality included in the 
# trajectory_analysis table.
        
### TABLE 17: ranking_global
            # Empty dictionary for the ranking
            Ranking = {}
            
            for file in glob.glob( PATH_RANKING + "SYSTEM" + "*" ):
                
                # Type of ranking
                kind = re.search( '_(.*)_', file.split("/")[-1] ).group(1)
                
                # Open the ranking file
                with open( file ) as FILE:
                    RANKING_LIST = json.load( FILE )
                    
                    # Find the position of the system in the ranking
                    for SIM in range(len(RANKING_LIST)):
                        if README["path"] in RANKING_LIST[SIM]["system"]["path"]: 
                            Ranking[ kind ] = SIM + 1
                    
                    # If it does not have an assigned value, use None
                    if not kind in Ranking:
                        Ranking[ kind ] = 4242
                
                FILE.close()
            
            # Collect the information of the position of the system in the ranking
            Info = { "trajectory_id": Trj_ID,
                     "ranking_total": Ranking["total"],
                     "ranking_hg":    Ranking["headgroup"],
                     "ranking_tails": Ranking["tails"],
                     "quality_total": QUALITY_SYSTEM["total"],
                     "quality_hg":    QUALITY_SYSTEM["headgroup"],
                     "quality_tails": QUALITY_SYSTEM["tails"] }
            
            # The minimal information about the system
            Minimal = { "trajectory_id": Trj_ID }
            
            # Entry in the DB with the info of ranking
            _ = DBEntry( 'ranking_global', Info, Minimal )
            
            
### TABLE 18: ranking_lipids
            # Empty dictionary for the ranking
            Ranking_lipids = {}
            
            for lipid in Lipids:
                Ranking_lipids[ lipid ] = {}
                
                for file in glob.glob( PATH_RANKING + lipid + "*" ):
                    # Type of ranking
                    kind = re.search( '_(.*)_', file ).group(1)
                    
                    # Open the ranking file
                    with open( file ) as FILE:
                        RANKING_LIST = json.load( FILE )
                        
                        # Find the position of the system in the ranking
                        for SIM in range(len(RANKING_LIST)):
                            if README["path"] in RANKING_LIST[SIM]["system"]["path"]: 
                                Ranking_lipids[ lipid ][ kind ] = SIM + 1
                    
                    FILE.close()
                    
                for t in ["total", "headgroup", "sn-1", "sn-2"]:
                    try:    Ranking_lipids[lipid][t] = Ranking_lipids[lipid][t] if not np.isnan( Ranking_lipids[lipid][t] ) else 4242
                    except: Ranking_lipids[lipid][t] = 4242
                        
                    #try:    Lipid_Quality[lipid][t] = Lipid_Quality[lipid][t]  if not np.isnan( Lipid_Quality[lipid][t] ) else 0
                    #except: Lipid_Quality[lipid][t] = 0
                
                # Collect the information of the position of the system in the ranking
                Info = { "trajectory_id": Trj_ID,
                         "lipid_id":      Lipids_ID[lipid],
                         "ranking_total": Ranking_lipids[lipid]["total"],
                         "ranking_hg":    Ranking_lipids[lipid]["headgroup"],
                         "ranking_sn-1":  Ranking_lipids[lipid]["sn-1"],
                         "ranking_sn-2":  Ranking_lipids[lipid]["sn-2"],
                         "quality_total": Lipid_Quality[lipid]["total"],
                         "quality_hg":    Lipid_Quality[lipid]["headgroup"],
                         "quality_sn-1":  Lipid_Quality[lipid]["sn-1"],
                         "quality_sn-2":  Lipid_Quality[lipid]["sn-2"] }
                
                # The minimal information about the system
                Minimal = { "trajectory_id": Trj_ID,
                            "lipid_id":      Lipids_ID[lipid] }
                
                # Entry in the DB with the info of ranking
                _ = DBEntry( 'ranking_lipids', Info, Minimal )
        
        
### TABLE 19: ranking_heteromolecules
            # Empty dictionary for the ranking
            Ranking_heteromolecules = {}
            
            for hetero in Heteromolecules:
                Ranking_heteromolecules[ hetero ] = {}
                
                for file in glob.glob( PATH_RANKING + hetero + "*" ):
                    # Type of ranking
                    kind = re.search( '_(.*)_', file ).group(1)
                    
                    # Open the ranking file
                    with open( file ) as FILE:
                        RANKING_LIST = json.load( FILE )
                        
                        # Find the position of the system in the ranking
                        for SIM in range(len(RANKING_LIST)):
                            if README["path"] in RANKING_LIST[SIM]["system"]["path"]: 
                                Ranking_heteromolecules[ hetero ][ kind ] = SIM + 1
                    
                    FILE.close()
                
                for t in ["total","headgroup","tail"]:
                    try:    Ranking_heteromolecules[hetero][t] = Ranking_heteromolecules[hetero][t] if not np.isnan( Ranking_heteromolecules[hetero][t] ) else 4242
                    except: Ranking_heteromolecules[hetero][t] = 4242
                        
                    #try:    Heteromolecules_Quality[hetero][t] = Heteromolecules_Quality[hetero][t] if not np.isnan( Heteromolecules_Quality[hetero][t] ) else 0
                    #except: Heteromolecules_Quality[hetero][t] = 0
                
                # Collect the information of the position of the system in the ranking
                Info = { "trajectory_id": Trj_ID,
                         "molecule_id":   Heteromolecules_ID[hetero],
                         "ranking_total": Ranking_heteromolecules[hetero]["total"],
                         "quality_total": Heteromolecules_Quality[hetero]["total"] }
                
                # The minimal information about the system
                Minimal = { "trajectory_id": Trj_ID,
                            "molecule_id":      Heteromolecules_ID[hetero] }
                
                # Entry in the DB with the info of ranking
                _ = DBEntry( 'ranking_heteromolecules', Info, Minimal )
                
                
        except:
            FAILS.append( README["path"] )

    if FAILS:
        print("\nThe following systems failed. Please check the files.")
        print( "\n" + "\n".join( FAILS ) )
        
####################

    database.close()
