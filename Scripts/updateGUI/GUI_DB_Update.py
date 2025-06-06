#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 14:27:05 2021

@author: Fabs

"""

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODULES
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import os
import os.path as osp
import re
import sys
import glob
import json
import yaml
import pymysql
import argparse
import numpy as np
import numbers
from importlib import import_module

# IMPORTLIB imports just `core` and `databankio` to avoid additional dependecies.
# It DOES NOT require the package to be pre-installed
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
dbl = import_module("DatabankLib")
NMRDict = import_module("DatabankLib.settings.molecules")
core = import_module("DatabankLib.core")
sys.path.pop(0)


# most of paths should be inserted into the DB relative to repo root
def genRpath(apath):
    return osp.relpath(apath, dbl.NMLDB_ROOT_PATH)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ARGUMENTS
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


# Program description
parser = argparse.ArgumentParser(description='NMRLipids Update v1.0')

# Ubication of data
parser.add_argument(
    "-c", "--config", type=str, default="config.json",
    help=''' JSON file with the configuration of the connection to the DB.
    Default: %(default)s ''')

# System properties
parser.add_argument(
    "-s", "--systems", type=str, nargs='+',  # REQUIRED
    help=""" Path of the system(s). """)

args = parser.parse_args()


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SQL Queries
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


def SQL_Select(Table: str, Values: list, Condition: dict = {}) -> str:
    '''
    Generate a SQL query to select values in a table. It compares floats with 1E-5
    tolerance!

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
          WHEN Condition.keys()[0]=Condition.value()[0] AND ...
               Condition.keys()[-1]=Condition.value()[-1]
    '''
    Query = (
        ' SELECT ' + ", ".join(map(lambda x: f'`{x}`', Values)) +
        f' FROM `{Table}` '
    )
    # Add a condition to the search
    if Condition:
        comps = []
        for k, v in Condition.items():
            if isinstance(v, numbers.Number) and v != np.ceil(v):
                comp = f'ABS( `{k}` - {v}) < 1E-5'
            else:
                comp = f'`{k}` = "{v}"'
            comps.append(comp)
        Query += 'WHERE ' + " AND ".join(comps)
    return Query


def SQL_Create(Table: str, Values: dict, Condition: dict = {}) -> str:
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
        INSERT INTO Table ( Values.keys()[0], ..., Values.keys()[-1] ) VALUES
                    ( Values.values()[0], ..., Values.values()[-1] )
               WHEN Condition.keys()[0]=Condition.value()[0] AND ...
                    Condition.keys()[-1]=Condition.value()[-1]
    '''
    Query = (
        f' INSERT INTO `{Table}` (' +
        ", ".join(map(lambda x: f'`{x}`', Values.keys())) +
        ") VALUES (" + ", ".join(map(lambda x: f'"{x}"', Values.values())) + ') '
    )

    # Add a condition to the search
    if Condition:
        Query += (
            'WHERE ' +
            " AND ".join(map(lambda x: f'`{x[0]}`="{x[1]}"', Condition.items()))
        )

    return Query


def SQL_Update(Table: str, Values: dict, Condition: dict = {}) -> str:
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
        UPDATE Table SET Values.keys()[0] = Values.values()[0], ...,
                         Values.keys()[-1] = Values.values()[-1]
          WHEN Condition.keys()[0]=Condition.value()[0] AND ...
               Condition.keys()[-1]=Condition.value()[-1]
    '''
    Query = (
        f' UPDATE `{Table}` SET ' +
        ', '.join(map(lambda x: f'`{x[0]}`="{x[1]}"', Values.items())) + ' '
    )

    # Add a condition to the search
    if Condition:
        Query += (
            'WHERE ' +
            " AND ".join(map(lambda x: f'`{x[0]}`="{x[1]}"', Condition.items()))
        )

    return Query


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Functions
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def CheckEntry(Table: str, Information: dict = {}) -> int:
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
    with database.cursor() as cursor:
        # Find the ID(s) of the entry matching the condition
        cursor.execute(SQL_Select(Table, ["id"], Information))
        # The list of IDs
        ID = cursor.fetchone()

    # return None imidiately if the record is not found
    if ID is None:
        return ID

    # More than ID will raise an error
    assert len(ID) == 1

    # extract ID-s value
    return ID[0]


def CreateEntry(Table: str, Information: dict) -> int:
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

    # Create a cursor
    with database.cursor() as cursor:
        # Execute the query creating a new entry
        res = cursor.execute(SQL_Create(Table, Information))

    # Commit the changes
    database.commit()

    # Num rows affected should be 1
    if res != 1:
        print("ERROR: record wasn't inserted!")
        print(Information)
        return 0

    # Check if the entry was created
    ID = CheckEntry(Table, Information)
    # If there is not an ID, raise an error (the table was not created)
    if not ID:
        print("WARNING: Something may have gone wrong with the table {}".format(Table))
        print(Information)
        return 0
    # If an ID is obtained, the entry was created succesfuly
    else:
        print("A new entry was created in {}: index {}".format(Table, ID))
        return ID


def UpdateEntry(Table: str, Information: dict, Condition: dict):
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
    with database.cursor() as cursor:
        # Execute the query updating an entry
        cursor.execute(SQL_Update(Table, Information, Condition))

    # Commit the changes
    database.commit()

    return print("Entry {} in table {} was updated".format(Condition["id"], Table))


def DBEntry(Table: str, Information: dict, Minimal: dict = {}) -> tuple:
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

    # --- TEMPORAL -----
    # Delete the Nan from the Information dictionary
    # The presence of NaN in the DB leads to problems dealing with the data
    # A solution must be found for this problem in the future.
    for entry in Information:
        if Information[entry] == "nan":
            Information[entry] = 0
    # --------------------

    # Check the existence of the entry
    EntryID = CheckEntry(Table, Minimal)

    # If there is an ID associated to the minimal information...
    if EntryID:

        # Check if the Minimal information matches the whole Information
        # If they don't match...
        if Information != Minimal:

            # Add the ID of the entry to the minimal information
            Minimal["id"] = EntryID

            # Update the entry
            UpdateEntry(Table, Information, Minimal)
            return EntryID

        # If the minimal information and the total information match,
        # no update is required.
        else:
            print("It was not necessary to update entry {} in table {}"
                  .format(EntryID, Table))
            return EntryID

    # If there is not an entry...
    else:
        # Create a new one with the data from Information
        EntryID = CreateEntry(Table, Information)
        return EntryID


def ResetTable(Table: str):
    '''
    Remove the content of a table and reset the index

    Parameters
    ----------
    Table : str
        Name of the table.

    Returns
    -------
    None.

    '''

    # Create a cursor
    with database.cursor() as cursor:
        cursor.execute(f' DELETE FROM `{Table}`')
        cursor.execute(f'ALTER TABLE `{Table}` AUTO_INCREMENT=1')

    return


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MAIN PROGRAM
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# The list of molecules in the membrane whose structure is not a phospholipid
HETEROMOLECULES_LIST = ["CHOL", "DCHOL", "C20", "C30"]
FAILS = []

TABLE_LIST = ["experiments_OP",
              "experiments_FF",
              "forcefields",
              "lipids",
              "ions",
              "water_models",
              "heteromolecules",
              "membranes",
              "trajectories",
              "trajectories_lipids",
              "trajectories_water",
              "trajectories_ions",
              "trajectories_heteromolecules",
              "trajectories_membranes",
              "trajectories_analysis",
              "trajectories_analysis_lipids",
              "trajectories_analysis_heteromolecules",
              "trajectories_analysis_ions",
              "trajectories_analysis_water",
              "ranking_global",
              "ranking_lipids",
              "ranking_heteromolecules",
              "trajectories_experiments_OP",
              "trajectories_experiments_FF"
              ]

if __name__ == '__main__':
    # Load the configuration of the connection
    config = json.load(open(args.config, "r"))

    # Conect to the database
    database = pymysql.connect(**config)

    # Path to the ranking files
    PATH_RANKING = osp.join(dbl.NMLDB_DATA_PATH, "Ranking")
    PATH_SIMULATIONS = dbl.NMLDB_SIMU_PATH
    PATH_EXPERIMENTS = dbl.NMLDB_EXP_PATH

    # If -s(ystems) is given, only modify that entries
    if args.systems is not None:
        systems = []
        for system in args.systems:
            # Format the system info
            system = system.strip("/\\")

            # The path to the simulation file
            PATH_SIMULATION = osp.join(PATH_SIMULATIONS, system)

            # Read the README.yaml file, the information of the simulation
            with open(osp.join(PATH_SIMULATION, 'README.yaml')) as File:
                README = yaml.load(File, Loader=yaml.FullLoader)

            # Add the math of the file
            README["path"] = system

            systems.append(README)

    else:
        # Initialize the whole databank
        systems = core.initialize_databank()

        # Remove the content of the tables and reset the autoincrements
        for table in TABLE_LIST[::-1]:
            ResetTable(table)

        # Create dummy FF and membrane
        Info = {
            "name":   'None',
            "date":   'None',
            "source": 'None'
            }
        _ = DBEntry('forcefields', Info, Info)

        # Create dummy FF and membrane
        Info = {
            "forcefield_id":   1,
            "lipid_names_l1":  "None",
            "lipid_names_l2":  "None",
            "lipid_number_l1": "None",
            "lipid_number_l2": "None",
            "geometry":        "None"
            }
        _ = DBEntry('membranes', Info, Info)


# -- TABLE `experiments_OP`

    # Find files with order parameters experiments
    EXP_OP = []
    PATH_EXPERIMENTS_OP = osp.join(PATH_EXPERIMENTS, "OrderParameters")

    # Get the path to every README.yaml file with experimental data
    for path, _, files in os.walk(PATH_EXPERIMENTS_OP):
        for file in files:
            if file.endswith("README.yaml"):
                EXP_OP.append(osp.relpath(path, PATH_EXPERIMENTS_OP))

    # Iterate over each experiment
    for expOP in EXP_OP:

        # Get the DOI of the experiment and the path to the README.yaml file
        with open(osp.join(PATH_EXPERIMENTS_OP, expOP, 'README.yaml')) as File:
            README = yaml.load(File, Loader=yaml.FullLoader)

        for file in os.listdir(osp.join(PATH_EXPERIMENTS_OP, expOP)):
            if file.endswith(".json"):

                Info = {"doi": README["DOI"],
                        "path": genRpath(osp.join(PATH_EXPERIMENTS_OP, expOP, file))}

                # Entry in the DB with the info of the experiment
                Exp_ID = DBEntry('experiments_OP', Info, Info)


# -- TABLE `experiments_FF`

    # Find files with form factor experiments
    EXP_FF = []
    PATH_EXPERIMENTS_FF = osp.join(PATH_EXPERIMENTS, "FormFactors")

    # Get the path to every README.yaml file with experimental data
    for path, _, files in os.walk(PATH_EXPERIMENTS_FF):
        for file in files:
            if file.endswith("README.yaml"):
                EXP_FF.append(osp.relpath(path, PATH_EXPERIMENTS_FF))

    # Iterate over each experiment
    for expFF in EXP_FF:

        # Get the DOI of the experiment and the path to the README.yaml file
        with open(osp.join(PATH_EXPERIMENTS_FF, expFF, 'README.yaml')) as File:
            README = yaml.load(File, Loader=yaml.FullLoader)

        for file in os.listdir(osp.join(PATH_EXPERIMENTS_FF, expFF)):
            if file.endswith(".json"):

                Info = {"doi": README["DOI"],
                        "path": genRpath(osp.join(PATH_EXPERIMENTS_FF, expFF, file))}

                # Entry in the DB with the info of the experiment
                Exp_ID = DBEntry('experiments_FF', Info, Info)

    # Iterate over the loaded systems
    for README in systems:

        try:
            # if True:
            print("\nCollecting data from system:")
            print(README["path"] + "\n")

            # The location of the files
            PATH_SIMULATION = osp.join(PATH_SIMULATIONS, README["path"])

            # In the case a field in the README does not exist, set its value to 0
            for field in [
                    'AUTHORS_CONTACT', 'COMPOSITION', 'CPT', 'DATEOFRUNNING', 'DIR_WRK',
                    'DOI', 'FF', 'FF_DATE', 'FF_SOURCE', 'GRO', 'LOG',
                    'NUMBER_OF_ATOMS', 'PREEQTIME', 'PUBLICATION', 'SOFTWARE',
                    'SOFTWARE_VERSION', 'SYSTEM', 'TEMPERATURE', 'TIMELEFTOUT', 'TOP',
                    'TPR', 'TRAJECTORY_SIZE', 'TRJ', 'TRJLENGTH', 'TYPEOFSYSTEM',
                    'WARNINGS', 'ID']:
                if field not in README:
                    README[field] = 0

    # -- TABLE `forcefields`
            # Collect the information about the forcefield
            Info = {
                "name":   README["FF"],
                "date":   README["FF_DATE"],
                "source": README["FF_SOURCE"]
                }

            # Entry in the DB with the info of the FF
            FF_ID = DBEntry('forcefields', Info, Info)

    # -- TABLE `lipids`
            # Empty dictionaries for the info of the lipids
            Lipids = {}
            Lipids_ID = {}
            Lipid_Ranking = {}
            Lipid_Quality = {}
            # Find the lipids in the composition
            for key in README["COMPOSITION"]:
                if key in NMRDict.lipids_set and key not in HETEROMOLECULES_LIST:
                    # Save the quality of the lipid
                    Store = True

                    # Collect the info of the lipids
                    Info = {
                        "forcefield_id": FF_ID,
                        "molecule":      key,
                        "name":          README["COMPOSITION"][key]["NAME"],
                        "mapping":       README["COMPOSITION"][key]["MAPPING"]
                        }

                    # Entry in the DB with the info of the lipid
                    Lip_ID = DBEntry('lipids', Info, Info)

                    # Store information for further steps
                    Lipids[key] = README["COMPOSITION"][key]["COUNT"]
                    Lipids_ID[key] = Lip_ID

    # --- TEMPORAL -----
    # Must be chenged when the final structure is ready
    # If the LIPID_FragmentQuality.json file will be defined for every system
    # this part can be deleted and just read the quality (try at the end of this
    # part). At this moment the ranking is not necessary in the DB, the web
    # already provides the result sorted by quality.
    # Same for the heteromolecules.
                    Lipid_Ranking[key] = {}
                    # Find the position of the system in the ranking
                    for file in glob.glob(osp.join(PATH_RANKING, key) + "*"):

                        # print(file)

                        # Kind of ranking (total, headgroup...)
                        kind = re.search('_(.*)_', file).group(1)

                        # Open the ranking file
                        with open(file) as FILE:
                            RANKING_LIST = json.load(FILE)

                            # Find the position of the system in the ranking
                            for SIM in range(len(RANKING_LIST)):

                                if README["path"] in \
                                        RANKING_LIST[SIM]["system"]["path"]:
                                    Lipid_Ranking[key][kind] = SIM + 1

                                    if Store:
                                        Lipid_Quality[key] = RANKING_LIST[SIM][key]
                                        Store = False

                            # If it does not have an assigned value, use None
                            if kind not in Lipid_Ranking[key]:
                                Lipid_Ranking[key][kind] = 0

                    # Read the quality file of the lipid (this will remain if the rest
                    # is removed)
                    try:
                        with open(osp.join(PATH_SIMULATION, key +
                                           '_FragmentQuality.json')) as FILE:
                            Lipid_Quality[key] = json.load(FILE)
                    except Exception:
                        Lipid_Quality[key] = {
                            "total": 0,
                            "headgroup": 0,
                            "sn-1": 0,
                            "sn-2": 0}

                    for t in ["total", "headgroup", "sn-1", "sn-2"]:
                        try:
                            Lipid_Quality[key][t] = Lipid_Quality[key][t]\
                                if not np.isnan(Lipid_Quality[key][t]) else 0
                        except Exception:
                            Lipid_Quality[key][t] = 0
    # ------------------

    # -- TABLE `ions`
            # Empty dictionary for the info of the ions
            Ions = {}
            # Find the ions in the composition
            for key in README["COMPOSITION"]:
                if key in NMRDict.molecules_set and key != "SOL" and \
                        key not in HETEROMOLECULES_LIST:
                    # Collect the info of the ions
                    Info = {
                        "forcefield_id": FF_ID,
                        "molecule":      key,
                        "name":          README["COMPOSITION"][key]["NAME"],
                        "mapping":       README["COMPOSITION"][key]["MAPPING"]
                        }

                    # Entry in the DB with the info of the ion
                    Ion_ID = DBEntry('ions', Info, Info)

                    # Store information for further steps: Ions[name]=[ID,number]
                    Ions[key] = [Ion_ID, README["COMPOSITION"][key]["COUNT"]]

    # -- TABLE `water_models`
            if "SOL" in README["COMPOSITION"]:
                WatName = README["COMPOSITION"]["SOL"]["NAME"]

                WatNum = README["COMPOSITION"]["SOL"]["COUNT"]

                # Collect the info on the water
                Info = {
                    "short_name":    WatName,
                    "mapping":       README["COMPOSITION"]["SOL"]["MAPPING"]
                    }

                # Get an entry in the DB with the info of the water
                Wat_ID = DBEntry('water_models', Info, Info)

    # -- TABLE `heteromolecules`
    # Heteromolecules are defined as the lipids for whom a distinction between the
    # different parts (headgroup, sn-1, sn-2) was not made. They could be included
    # in the lipids table and leaving some fields empty.

            # Empty dictionaries for the info of the heteromolecules
            Heteromolecules = {}
            Heteromolecules_ID = {}
            Heteromolecules_Ranking = {}

            # Form factor
            Store = True
            Heteromolecules_Quality = {}

            # Find the heteromolecules in the composition
            for key in README["COMPOSITION"]:
                if key in NMRDict.lipids_set and key in HETEROMOLECULES_LIST:

                    # Collect the info of the lipids
                    Info = {
                        "forcefield_id": FF_ID,
                        "molecule":      key,
                        "name":          README["COMPOSITION"][key]["NAME"],
                        "mapping":       README["COMPOSITION"][key]["MAPPING"]
                        }

                    # Entry in the DB with the info of the heteromolecules
                    Mol_ID = DBEntry('heteromolecules', Info, Info)

                    # Store information for further steps
                    Heteromolecules[key] = README["COMPOSITION"][key]["COUNT"]
                    Heteromolecules_ID[key] = Mol_ID

    # --- TEMPORAL -----
    # See the lipids table for the reasons
                    Heteromolecules_Ranking[key] = {}

                    # Find the position of the system in the raking
                    for file in glob.glob(osp.join(PATH_RANKING, key) + "*"):

                        # Type of ranking
                        kind = re.search('_(.*)_', file).group(1)

                        # Open the ranking file
                        with open(file) as FILE:
                            RANKING_LIST = json.load(FILE)

                            # Find the position of the system in the ranking
                            for SIM in range(len(RANKING_LIST)):
                                if README["path"] in \
                                        RANKING_LIST[SIM]["system"]["path"]:
                                    Heteromolecules_Ranking[key][kind] = SIM + 1

                                    if Store:
                                        Heteromolecules_Quality[key] = \
                                            RANKING_LIST[SIM][key]
                                        Store = False

                            # If it does not have an assigned value, use None
                            if kind not in Heteromolecules_Ranking[key]:
                                Heteromolecules_Ranking[key][kind] = 0

                    try:
                        with open(osp.join(PATH_SIMULATION,
                                           key + '_FragmentQuality.json')) as FILE:
                            Heteromolecules_Quality[key] = json.load(FILE)
                    except Exception:
                        Heteromolecules_Quality[key] = {
                            "total": 0,
                            "headgroup": 0,
                            "tail": 0
                            }

                    for t in ["total", "headgroup", "tail"]:
                        try:
                            Heteromolecules_Quality[key][t] = \
                                Heteromolecules_Quality[key][t] \
                                if not np.isnan(Heteromolecules_Quality[key][t])\
                                else 0
                        except Exception:
                            Heteromolecules_Quality[key][t] = 0
    # ------------------

    # -- TABLE `membranes`
            # Find the proportion of each lipid in the leaflets
            Names = [[], []]
            Number = [[], []]

            for lipid in Lipids:
                if Lipids[lipid][0]:
                    Names[0].append(lipid)
                    Number[0].append(str(Lipids[lipid][0]))
                if Lipids[lipid][1]:
                    Names[1].append(lipid)
                    Number[1].append(str(Lipids[lipid][1]))

            for hetero in Heteromolecules:
                if Heteromolecules[hetero][0]:
                    Names[0].append(hetero)
                    Number[0].append(str(Heteromolecules[hetero][0]))
                if Heteromolecules[hetero][1]:
                    Names[1].append(hetero)
                    Number[1].append(str(Heteromolecules[hetero][1]))

            Names = [':'.join(Names[0]), ':'.join(Names[1])]
            Number = [':'.join(Number[0]), ':'.join(Number[1])]

            # Collect the information about the membrane
            Info = {
                "forcefield_id":   FF_ID,
                "lipid_names_l1":  Names[0],
                "lipid_names_l2":  Names[1],
                "lipid_number_l1": Number[0],
                "lipid_number_l2": Number[1],
                "geometry":        README["TYPEOFSYSTEM"]
                }

            # Entry in the DB with the info of the membrane
            Mem_ID = DBEntry('membranes', Info, Info)

    # -- TABLE `trajectories`
            # Collect the information about the simulation
            Info = {
                "id":              README["ID"],
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
                "trj_length":      README["TRJLENGTH"]
                }

            # The information that defines the trajectory
            Minimal = {
                "id":            README["ID"],
                "forcefield_id": FF_ID,
                "membrane_id":   Mem_ID,
                "git_path":      README["path"],
                "system":        README["SYSTEM"]
                }

            # Entry in the DB with the info of the trajectory
            Trj_ID = DBEntry('trajectories', Info, Minimal)

    # -- TABLE `trajectories_lipids`
            TrjL_ID = {}
            for lipid in Lipids:
                # Collect the information of each lipid in the simulation
                Info = {
                    "trajectory_id": Trj_ID,
                    "lipid_id":      Lipids_ID[lipid],
                    "lipid_name":    lipid,
                    "leaflet_1":     Lipids[lipid][0],
                    "leaflet_2":     Lipids[lipid][1]
                    }

                # The minimal information that identifies the lipid
                Minimal = {
                    "trajectory_id": Trj_ID,
                    "lipid_id":      Lipids_ID[lipid]
                    }

                # Entry in the DB with the info of the lipids in the simulation
                TrjL_ID[lipid] = DBEntry('trajectories_lipids', Info, Minimal)

    # -- TABLE `trajectories_water`
            if "SOL" in README["COMPOSITION"]:
                # Collect the information of the water in the simulation
                Info = {
                    "trajectory_id": Trj_ID,
                    "water_id":      Wat_ID,
                    "water_name":    WatName,
                    "number":        WatNum}

                # The minimal information of the water in the simulation
                Minimal = {
                    "trajectory_id": Trj_ID,
                    "water_id":      Wat_ID}

                # Entry in the DB with the info of the water in the simulation
                _ = DBEntry('trajectories_water', Info, Minimal)

    # -- TABLE `trajectories_ions`
            TrjI_ID = {}
            for ion in Ions:
                # Collect the information of each ion in the simulation
                Info = {
                    "trajectory_id": Trj_ID,
                    "ion_id":        Ions[ion][0],
                    "ion_name":      ion,
                    "number":        Ions[ion][1]}

                # The minimal information that identifies the ion
                Minimal = {
                    "trajectory_id": Trj_ID,
                    "ion_id":        Ions[ion][0]}

                # Entry in the DB with the info of the ions in the simulation
                TrjI_ID[ion] = DBEntry('trajectories_ions', Info, Minimal)

    # -- TABLE `trajectories_heteromolecules`
            TrjM_ID = {}
            for hetero in Heteromolecules:
                # Collect the information of each heteromolecule in the simulation
                Info = {
                    "trajectory_id": Trj_ID,
                    "molecule_id":   Heteromolecules_ID[hetero],
                    "molecule_name": hetero,
                    "leaflet_1":     Heteromolecules[hetero][0],
                    "leaflet_2":     Heteromolecules[hetero][1]
                    }

                # The minimal information that identifies the heteromolecule
                Minimal = {
                    "trajectory_id": Trj_ID,
                    "molecule_id":   Heteromolecules_ID[hetero]}

                # Entry in the DB with the info of the heteromolecules in the simulation
                TrjM_ID[hetero] = DBEntry('trajectories_heteromolecules', Info, Minimal)

    # -- TABLE `trajectories_membranes``

            Info = {
                "trajectory_id": Trj_ID,
                "membrane_id": Mem_ID,
                "name": README["SYSTEM"]}

            _ = DBEntry('trajectories_membranes', Info, Info)

    # -- TABLE `trajectories_analysis`
            # Find the bilayer thickness
            try:
                with open(osp.join(PATH_SIMULATION, 'thickness.json')) as FILE:
                    BLT = json.load(FILE)
            except Exception:
                BLT = 0

            # Find the area per lipid
            try:
                with open(osp.join(PATH_SIMULATION, 'apl.json')) as FILE:
                    # Load the file
                    ApL = json.load(FILE)

                    # Transform the dictionary into an array
                    ApL = np.array([[float(key), float(ApL[key])] for key in ApL])

                    # Perform the mean
                    APL = np.mean(ApL[int(len(ApL[:, 0])/2):, 1])
            except Exception:
                APL = 0

            # Form factor quality
            try:
                with open(osp.join(PATH_SIMULATION, 'FormFactorQuality.json')) as FILE:
                    FFQ = json.load(FILE)
            except Exception:
                FFQ = [4242, 0]

            # Read the quality file for the whole system
            try:
                with open(osp.join(PATH_SIMULATION, 'SYSTEM_quality.json')) as FILE:
                    QUALITY_SYSTEM = json.load(FILE)
            except Exception:
                QUALITY_SYSTEM = {
                    "total": 0,
                    "headgroup": 0,
                    "tails": 0}

            try:
                FFExp = genRpath(
                    osp.join(PATH_EXPERIMENTS_FF, README["EXPERIMENT"]["FORMFACTOR"]))
            except Exception:
                FFExp = ''

            # Collect the information of the analysis of the trajectory
            Info = {
                "trajectory_id":          Trj_ID,
                "bilayer_thickness":      BLT,
                "area_per_lipid":         APL,
                "area_per_lipid_file":    genRpath(
                    osp.join(PATH_SIMULATIONS, README["path"], 'apl.json')),
                "form_factor_file":       genRpath(
                    osp.join(PATH_SIMULATIONS, README["path"], 'FormFactor.json')),
                "quality_total":          QUALITY_SYSTEM["total"],
                "quality_headgroups":     QUALITY_SYSTEM["headgroup"],
                "quality_tails":          QUALITY_SYSTEM["tails"],
                "form_factor_experiment": FFExp,
                "form_factor_quality":    FFQ[0],
                "form_factor_scaling":    FFQ[1]
                }

            # Collect the minimal information of the analysis of the trajectory
            Minimal = {"trajectory_id": Trj_ID}

            # Entry in the DB with the info of the analysis of the simulation
            _ = DBEntry('trajectories_analysis', Info, Minimal)

    # -- TABLE `trajectories_analysis_lipids`
            for lipid in Lipids:
                try:
                    OPExp = genRpath(osp.join(
                        PATH_EXPERIMENTS_OP,
                        list(README["EXPERIMENT"]["ORDERPARAMETER"][lipid].values())[0],
                        lipid + '_Order_Parameters.json')
                        )
                except Exception:
                    OPExp = ''

                # Collect the information of each lipid in the simulation
                Info = {
                    "trajectory_id":                Trj_ID,
                    "lipid_id":                     Lipids_ID[lipid],
                    "quality_total":                Lipid_Quality[lipid]["total"],
                    "quality_hg":                   Lipid_Quality[lipid]["headgroup"],
                    "quality_sn-1":                 Lipid_Quality[lipid]["sn-1"],
                    "quality_sn-2":                 Lipid_Quality[lipid]["sn-1"],
                    "order_parameters_file":        genRpath(
                        osp.join(PATH_SIMULATIONS, README["path"],
                                 lipid + 'OrderParameters.json')),
                    "order_parameters_experiment":  OPExp,
                    "order_parameters_quality":     genRpath(
                        osp.join(PATH_SIMULATIONS, README["path"],
                                 lipid + '_OrderParameters_quality.json'))
                    }

                # The minimal information that identifies the lipid in the simulation
                Minimal = {"trajectory_id": Trj_ID,
                           "lipid_id":      Lipids_ID[lipid]}

                # Entry in the DB with the info of the analysis of the lipid
                # in the simulation
                _ = DBEntry('trajectories_analysis_lipids', Info, Minimal)

    # -- TABLE `trajectories_analysis_heteromolecules`
            for hetero in Heteromolecules:
                try:
                    OPExp = genRpath(
                        osp.join(
                            PATH_EXPERIMENTS_OP,
                            list(README["EXPERIMENT"]["ORDERPARAMETER"][hetero]
                                 .values())[0],
                            hetero + '_Order_Parameters.json')
                        )
                except Exception:
                    OPExp = ''

                # Collect the information of each heteromolecule in the simulation
                Info = {
                    "trajectory_id":   Trj_ID,
                    "molecule_id":     Heteromolecules_ID[hetero],
                    "quality_total":   Heteromolecules_Quality[hetero]["total"],
                    "quality_hg":      Heteromolecules_Quality[hetero]["headgroup"],
                    "quality_tails":   Heteromolecules_Quality[hetero]["tail"],
                    "order_parameters_file": genRpath(
                        osp.join(PATH_SIMULATIONS, README["path"],
                                 hetero + 'OrderParameters.json')),
                    "order_parameters_experiment":  OPExp,
                    "order_parameters_quality":     genRpath(
                        osp.join(PATH_SIMULATIONS, README["path"],
                                 hetero + '_OrderParameters_quality.json'))
                    }

                # The minimal information that identifies the heteromolecule in the
                # simulation
                Minimal = {"trajectory_id": Trj_ID,
                           "molecule_id":   Heteromolecules_ID[hetero]}

                # Entry in the DB with the info of the analysis of the heteromolecule
                # in the simulation
                _ = DBEntry('trajectories_analysis_heteromolecules', Info, Minimal)

    # -- TABLE `trajectory_analysis_ions`
            for ion in Ions:
                # Collect the information of the ions in the simulation
                Info = {"trajectory_id": Trj_ID,
                        "ion_id":        Ions[ion][0]}

                # The minimal information that identifies the ion in the simulation
                # Minimal = { "trajectory_id": Trj_ID,
                #            "ion_id":        Ions[ ion ][0] }

                # Entry in the DB with the info of the analysis of the ion in the
                # simulation
                _ = DBEntry('trajectories_analysis_ions', Info, Info)

    # -- TABLE `trajectory_analysis_water`
            if "SOL" in README["COMPOSITION"]:
                # Collect the information of the water in the simulation
                Info = {"trajectory_id": Trj_ID,
                        "water_id":      Wat_ID}
                # -     "density_file":  args.densities_folder + README["path"]
                # -                       + '/SOLdensity.xvg' }

                # Entry in the DB with the info of the analysis of the water in the
                # simulation
                _ = DBEntry('trajectories_analysis_water', Info, Info)

    # --- TEMPORAL -----
    # The table may be removed in the future, and the quality included in the
    # trajectory_analysis table.

    # -- TABLE `ranking_global`
            # Empty dictionary for the ranking
            Ranking = {}

            for file in glob.glob(osp.join(PATH_RANKING, "SYSTEM") + "*"):

                # Type of ranking
                kind = re.search('_(.*)_', file.split("/")[-1]).group(1)

                # Open the ranking file
                with open(file) as FILE:
                    RANKING_LIST = json.load(FILE)

                    # Find the position of the system in the ranking
                    for SIM in range(len(RANKING_LIST)):
                        if README["path"] in RANKING_LIST[SIM]["system"]["path"]:
                            Ranking[kind] = SIM + 1

                    # If it does not have an assigned value, use None
                    if kind not in Ranking:
                        Ranking[kind] = 4242

            # Collect the information of the position of the system in the ranking
            Info = {
                "trajectory_id": Trj_ID,
                "ranking_total": Ranking["total"],
                "ranking_hg":    Ranking["headgroup"],
                "ranking_tails": Ranking["tails"],
                "quality_total": QUALITY_SYSTEM["total"],
                "quality_hg":    QUALITY_SYSTEM["headgroup"],
                "quality_tails": QUALITY_SYSTEM["tails"]
                }

            # The minimal information about the system
            Minimal = {"trajectory_id": Trj_ID}

            # Entry in the DB with the info of ranking
            _ = DBEntry('ranking_global', Info, Minimal)

    # -- TABLE `ranking_lipids`
            # Empty dictionary for the ranking
            Ranking_lipids = {}

            for lipid in Lipids:
                Ranking_lipids[lipid] = {}

                for file in glob.glob(osp.join(PATH_RANKING, lipid) + "*"):
                    # Type of ranking
                    kind = re.search('_(.*)_', file).group(1)

                    # Open the ranking file
                    with open(file) as FILE:
                        RANKING_LIST = json.load(FILE)

                        # Find the position of the system in the ranking
                        for SIM in range(len(RANKING_LIST)):
                            if README["path"] in RANKING_LIST[SIM]["system"]["path"]:
                                Ranking_lipids[lipid][kind] = SIM + 1

                for t in ["total", "headgroup", "sn-1", "sn-2"]:
                    try:
                        Ranking_lipids[lipid][t] = Ranking_lipids[lipid][t] \
                            if not np.isnan(Ranking_lipids[lipid][t]) else 4242
                    except Exception:
                        Ranking_lipids[lipid][t] = 4242

                # Collect the information of the position of the system in the ranking
                Info = {
                    "trajectory_id": Trj_ID,
                    "lipid_id":      Lipids_ID[lipid],
                    "ranking_total": Ranking_lipids[lipid]["total"],
                    "ranking_hg":    Ranking_lipids[lipid]["headgroup"],
                    "ranking_sn-1":  Ranking_lipids[lipid]["sn-1"],
                    "ranking_sn-2":  Ranking_lipids[lipid]["sn-2"],
                    "quality_total": Lipid_Quality[lipid]["total"],
                    "quality_hg":    Lipid_Quality[lipid]["headgroup"],
                    "quality_sn-1":  Lipid_Quality[lipid]["sn-1"],
                    "quality_sn-2":  Lipid_Quality[lipid]["sn-2"]
                    }

                # The minimal information about the system
                Minimal = {"trajectory_id": Trj_ID,
                           "lipid_id":      Lipids_ID[lipid]}

                # Entry in the DB with the info of ranking
                _ = DBEntry('ranking_lipids', Info, Minimal)

    # -- TABLE `ranking_heteromolecules`
            # Empty dictionary for the ranking
            Ranking_heteromolecules = {}

            for hetero in Heteromolecules:
                Ranking_heteromolecules[hetero] = {}

                for file in glob.glob(osp.join(PATH_RANKING, hetero) + "*"):
                    # Type of ranking
                    kind = re.search('_(.*)_', file).group(1)

                    # Open the ranking file
                    with open(file) as FILE:
                        RANKING_LIST = json.load(FILE)

                        # Find the position of the system in the ranking
                        for SIM in range(len(RANKING_LIST)):
                            if README["path"] in RANKING_LIST[SIM]["system"]["path"]:
                                Ranking_heteromolecules[hetero][kind] = SIM + 1

                for t in ["total", "headgroup", "tail"]:
                    try:
                        Ranking_heteromolecules[hetero][t] = \
                            Ranking_heteromolecules[hetero][t] \
                            if not np.isnan(Ranking_heteromolecules[hetero][t])\
                            else 4242
                    except Exception:
                        Ranking_heteromolecules[hetero][t] = 4242

                # Collect the information of the position of the system in the ranking
                Info = {
                    "trajectory_id": Trj_ID,
                    "molecule_id":   Heteromolecules_ID[hetero],
                    "ranking_total": Ranking_heteromolecules[hetero]["total"],
                    "quality_total": Heteromolecules_Quality[hetero]["total"]
                    }

                # The minimal information about the system
                Minimal = {"trajectory_id": Trj_ID,
                           "molecule_id":      Heteromolecules_ID[hetero]}

                # Entry in the DB with the info of ranking
                _ = DBEntry('ranking_heteromolecules', Info, Minimal)

            if "EXPERIMENT" in README:
                if "ORDERPARAMETER" in README["EXPERIMENT"]:
                    # -- TABLE `trajectories_experiments_OP`
                    # The Order Parameters experiments associated to the simulation

                    ExpOP = README["EXPERIMENT"]["ORDERPARAMETER"]
                    # Iterate over the lipids
                    for mol in ExpOP:
                        # Check if there is an experiment associated to the lipid
                        if ExpOP[mol]:
                            for doi, path in ExpOP[mol].items():
                                for file in os.listdir(
                                        osp.join(PATH_EXPERIMENTS_OP, path)):
                                    if file.endswith(".json"):
                                        Info = {
                                            "trajectory_id": Trj_ID,
                                            "lipid_id": {**Lipids_ID,
                                                         **Heteromolecules_ID}[mol],
                                            "experiment_id":
                                                CheckEntry(
                                                  'experiments_OP', {
                                                      "doi": doi,
                                                      "path": genRpath(osp.join(
                                                          PATH_EXPERIMENTS_OP,
                                                          path, file))
                                                        }
                                                    )
                                               }

                                        _ = DBEntry('trajectories_experiments_OP',
                                                    Info, Info)

    # -- TABLE `trajectories_experiments_FF`
                if "FORMFACTOR" in README["EXPERIMENT"]:
                    # The Form Factor experiments associated to the simulation
                    ExpFF = README["EXPERIMENT"]["FORMFACTOR"]

                    if ExpFF:
                        if type(ExpFF) is str:
                            ExpFF = [ExpFF]

                            for path in ExpFF:

                                for file in os.listdir(osp.join(
                                        PATH_EXPERIMENTS_FF, path)):
                                    if file.endswith(".json"):
                                        Info = {
                                            "trajectory_id": Trj_ID,
                                            "experiment_id": CheckEntry(
                                                     'experiments_FF', {
                                                         "path": genRpath(osp.join(
                                                             PATH_EXPERIMENTS_FF,
                                                             path, file))
                                                     })
                                                 }

                                        _ = DBEntry('trajectories_experiments_FF',
                                                    Info, Info)

        except Exception as err:
            print("Exception: "+str(err))
            FAILS.append(README["path"])

    # -- TABLE `trajectories` (again)
    try:

        # Generate a new cursor
        with database.cursor() as cursor:

            # Get the list of IDs currently in the DB
            cursor.execute(SQL_Select("trajectories", ["id"]))
            List_IDs = cursor.fetchall()

            maxID = int(open(
                osp.join(PATH_SIMULATIONS, "COUNTER_ID"), "r").readlines()[0])

            missing_IDs = set(range(1, maxID+1)) - {ID[0] for ID in List_IDs}

            Info = {
                "id":              1,
                "forcefield_id":   1,
                "membrane_id":     1,
                "git_path":        "''",
                "system":          "''",
                "author":          "''",
                "date":            "''",
                "dir_wrk":         "''",
                "doi":             "''",
                "number_of_atoms": 0,
                "preeq_time":      "''",
                "publication":     "''",
                "temperature":     0,
                "software":        "''",
                "trj_size":        0,
                "trj_length":      0,
                "timeleftout":     0
                }

            for missing_ID in missing_IDs:

                print("\n Adding missing ID:", missing_ID)
                Info["id"] = missing_ID

                cursor.execute(
                    "INSERT INTO `trajectories` (" +
                    ','.join(map(lambda x: '`' + str(x) + '`', Info.keys())) +
                    ") VALUES (" +
                    ','.join(map(str, Info.values())) +
                    ") "
                    )

        database.commit()

    except Exception:
        pass

    if FAILS:
        print(
            "\nThe following systems failed. Please check the files." +
            "\n" + "\n".join(FAILS)
            )

####################

    database.close()
