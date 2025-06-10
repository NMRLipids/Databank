#!/usr/bin/env python3
# coding: utf-8

"""
:program: AddData.py
:description: The script adds a simulation into the Databank based on ``info.yaml``
              file.

Usage:
    AddData.py Script [-h] [-f FILE] [-d] [-n] [-w WORK_DIR] [-o OUTPUT_DIR]

:param: ``-h``, ``--help`` show this help message and exit
:param: ``-f FILE``, ``--file``  Input config file in yaml format.
:param: ``-d``, ``--debug`` enable debug logging output
:param: ``-n``, ``--no-cache`` always redownload repository files
:param: ``-w WORK_DIR``, ``--work-dir`` set custom temporary working directory
        [not set = read from YAML]
:param: ``-o OUTPUT_DIR``, ``--output-dir`` set custom output directory
        [/path-to-databank]

Returning error codes:
    1 - input YAML parsing errors
    2 - filesystem writting errors
    3 - network accessing errors
"""

import os
import argparse
import yaml
import logging
import shutil
import pprint
import traceback
from datetime import date
from random import randint
from urllib.error import URLError, HTTPError
from copy import deepcopy
import pandas as pd
import numpy as np

from MDAnalysis import Universe

# import databank dictionaries
import DatabankLib
from DatabankLib.core import System
from DatabankLib.databankLibrary import (
    calc_file_sha1_hash,
    create_databank_directories,
    lipids_set,
    molecules_set
)
# helpers
from DatabankLib.databankio import (
    download_resource_from_uri,
    resolve_download_file_url
)
from DatabankLib.databankLibrary import (
    parse_valid_config_settings
)
from DatabankLib.settings.molecules import Lipid, NonLipid
from DatabankLib.settings.engines import get_struc_top_traj_fnames, software_dict

pd.set_option("display.max_rows", 500)
pd.set_option("display.max_columns", 500)
pd.set_option("display.width", 1000)
pd.set_option("display.max_colwidth", 1000)


if __name__ == "__main__":
    # parse input yaml file
    parser = argparse.ArgumentParser(
        prog="AddData.py Script",
        description="Add a new dataset to the NMRLipids databank"
    )
    parser.add_argument(
        "-f", "--file",
        help="Input config file in yaml format."
    )
    parser.add_argument(
        "-d", "--debug",
        help="enable debug logging output",
        action="store_true"
    )
    parser.add_argument(
        "-n", "--no-cache",
        help="always redownload repository files",
        action="store_true"
    )
    parser.add_argument(
        "-w", "--work-dir",
        help="set custom temporary working directory [not set = read from YAML]",
        default=""
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        help=f"set custom output directory [{DatabankLib.NMLDB_SIMU_PATH}]",
        default=DatabankLib.NMLDB_SIMU_PATH
    )

    args = parser.parse_args()

    # configure logging
    logging_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s",
        datefmt="%I:%M:%S %p",
        level=logging_level,
    )
    logger = logging.getLogger()

    input_path = os.path.normpath(args.file)

    # load input yaml file into empty dictionary
    info_yaml = {}

    # open input file for reading and writing
    with open(input_path) as yaml_file:
        info_yaml = yaml.load(
            yaml_file, Loader=yaml.FullLoader
        )  # TODO may throw yaml.YAMLError

    # Show the input read
    logger.debug(f"{os.linesep} Input read from {input_path} file:")
    pp = pprint.PrettyPrinter(width=41, compact=True)
    if logger.isEnabledFor(logging.DEBUG):
        pp.pprint(yaml.dump(info_yaml))

    # validate yaml entries and return updated sim dict
    try:
        sim_dict, files = parse_valid_config_settings(info_yaml)
    except KeyError as e:
        logger.error(f"missing entry key in yaml config: {e}, aborting")
        logger.error(traceback.format_exc())
        quit(1)
    except Exception as e:
        logger.error(
            f"an '{type(e).__name__}' occured while performing validity check"
            f" '{input_path}', script has been aborted"
        )
        logger.error(e)
        quit(1)
    else:
        logger.info(
            "all entries in simulation are understood and will be further processed"
        )
        logger.debug("valid sim entry keys:")
        pp = pprint.PrettyPrinter(width=41, compact=True)
        if logger.isEnabledFor(logging.DEBUG):
            pp.pprint(sim_dict)

    try:
        sim = System(sim_dict)
        # mapping files are registered here!
    except Exception as e:
        logger.error(
            f"an '{type(e).__name__}' occured while processing dict->System"
            f" '{input_path}', script has been aborted"
        )
        logger.error(e)
        quit(1)
    else:
        logger.info(f"System object is successfully created from '{input_path}' file")

    # Create temporary directory where to download files and analyze them

    if args.work_dir:
        dir_wrk = args.work_dir
        logger.warning(
            f"--work_dir override, ignoring 'DIR_WRK' from "
            f"configuration file: {sim['DIR_WRK']}"
        )
    else:
        dir_wrk = sim["DIR_WRK"]

    dir_tmp = (
        os.path.join(dir_wrk, "tmp_6-" + str(randint(100000, 999999)))
        if args.no_cache
        else os.path.join(dir_wrk, f"{sim['DOI'].split('/')[-1]}_download")
    )

    logger.info(f"The data will be processed in directory path '{dir_tmp}'")

    try:
        os.makedirs(dir_tmp, exist_ok=True)
    except OSError as e:
        logger.error(
            f"couldn't create temporary working directory '{dir_tmp}': {e.args[1]}"
        )
        quit(2)

    # Check link status and download files

    try:
        download_links = []
        for fi in files:
            logger.info(f"Validating file: {fi}..")
            _x = resolve_download_file_url(sim["DOI"], fi, validate_uri=True)
            download_links.append(_x)

        logger.info(f"Now downloading {len(files)} files ...")

        for url, fi in zip(download_links, files):
            download_resource_from_uri(
                url, os.path.join(dir_tmp, fi), override_if_exists=args.no_cache
            )

        logger.info(f"Download of {len(files)} files was successful")

    except HTTPError as e:
        if e.code == 404:
            logger.error(
                f"ressource not found on server '{e.url}' (404)."
                " Wrong DOI link or file name?"
            )
        else:
            logger.error(f"HTTPError {e.code} while trying to "
                         f"download the file '{e.url}'")
        quit(3)
    except URLError as e:
        logger.error(
            f"couldn't resolve network adress: {e.reason}."
            " Please check your internet connection."
        )
        quit(3)
    except Exception as e:
        logger.error(
            f"'{type(e).__name__}' while attempting to download ressources, aborting"
        )
        logger.error(traceback.format_exc())
        quit(3)

    # -- Calculate hash of downloaded files

    sim_hashes = deepcopy(sim)

    software_sim = software_dict[sim["SOFTWARE"].upper()]

    # list_containing the sha1 sums for all required files
    sha1_list_requied = []

    # Make empty dataframe with the desired columns
    df_files = pd.DataFrame(
        columns=["NAME", "TYPE", "REQUIRED", "SIZE_MB", "HASH"], dtype=object
    )

    for key_sim, value_sim in sim_hashes.items():
        # double-checking keys
        try:
            entry_type = software_sim[key_sim]["TYPE"]
        except KeyError:
            if key_sim in ["SOFTWARE", "ID"]:
                continue
            else:
                # That shouldn't happen! Unexpected YAML-keys were checked by
                # parse_valid_config_settings before
                raise

        if "file" in entry_type:
            files_list = []
            is_required = software_dict[sim_hashes["SOFTWARE"].upper()][key_sim][
                "REQUIRED"
            ]

            if not is_required and value_sim is None:
                continue  # skip not required NoneType (empty) file entries

            for file_provided in value_sim:
                file_name = os.path.join(dir_tmp, file_provided[0])
                logger.info(f"calculating sha1 hash of '{file_provided[0]}'...")
                file_hash = calc_file_sha1_hash(file_name)
                file_size_mb = f"{(os.path.getsize(file_name)/1024/1024):.2f}"

                df_files = pd.concat(
                    [
                        df_files,
                        pd.DataFrame(
                            [
                                {
                                    "NAME": file_provided[0],
                                    "TYPE": key_sim,
                                    "REQUIRED": is_required,
                                    "SIZE_MB": file_size_mb,
                                    "HASH": file_hash,
                                }
                            ]
                        ),
                    ],
                    ignore_index=True,
                )
                files_list.append([file_provided[0], file_hash])

                # Find the keys of the required files to calculate the master_hash
                if is_required:
                    sha1_list_requied.append(file_hash)

                sim_hashes[key_sim] = files_list  # TODO Problematic

    logger.info(f"Summary of downloaded files: {os.linesep}")
    logger.info("\n" + df_files.to_string())

    # Calculate the hash of a file contaning the hashes of each of the required files
    # This should be always invariant as it will be used unique identifier for a
    # simualtion. Note order the hashes of the required files before calculating the
    # hash (That means that the required files cannot change)

    # Calculates numbers of lipid molecules in each leaflet. This is done by checking
    # on which side of the centre of mass the membrane each the centre of mass of a
    # lipid molecule is. If a lipid molecule is split so that headgroup and tails are
    # their own residues, the centre of mass of the headgroup is used in the
    # calculation.
    ####################################################################################

    logger.info(
        "Calculating the numbers of lipid molecules in each leaflet based on the "
        "center of mass of the membrane and lipids."
    )
    logger.info(
        "If a lipid molecule is split to multiple residues, the centre of mass of"
        " the headgroup is used."
    )

    try:
        struc, top, traj = get_struc_top_traj_fnames(sim, join_path=dir_tmp)
    except (ValueError, KeyError) as e:
        logger.error(str(type(e)) + " => " + str(e))
        quit(1)
    except Exception as e:
        logger.error("Unkonwn error during `get_3major_fnames..`")
        logger.error(str(type(e)) + " => " + str(e))
        quit(4)

    leaflet1 = 0  # total number of lipids in upper leaflet
    leaflet2 = 0  # total number of lipids in lower leaflet

    # try to generate zero-frame structure from top + trajectory
    fail_from_top = False
    gro = os.path.join(dir_tmp, "frame0.gro")  # structure regenerated from top
    try:
        logger.info(f"MDAnalysis tries to use {top} and {traj}")
        u = Universe(top, traj)
        u.atoms.write(gro, frames=u.trajectory[[0]])
    except Exception as e:
        logger.warning(str(type(e)) + " => " + str(e))
        fail_from_top = True

    # if previous fails then try the same from struc + trajectory
    if fail_from_top and struc is not None:
        try:
            logger.info(f"MDAnalysis tries to use {struc} and {traj}")
            u = Universe(struc, traj)
            u.atoms.write(gro, frames=u.trajectory[[0]])
        except Exception as e:
            logger.error("Cannot initialize MDAnalysis using given structure file!")
            logger.error(str(type(e)) + " => " + e)
            quit(2)

    # if there is no struc and MDAnalysis doesn't start from TOP, then
    # GROMACS can try making struc from top!
    if fail_from_top and struc is None and sim["SOFTWARE"].upper() == "GROMACS":
        logger.info(
            "Now generating frame0.gro with Gromacs because MDAnalysis cannot "
            "read tpr version ..."
        )
        if (
            "WARNINGS" in sim and
            sim["WARNINGS"] is not None and
            sim["WARNINGS"]["GROMACS_VERSION"] == "gromacs3"
           ):
            exec_str = (
                f"executing 'echo System | trjconv -s {top} -f {traj} "
                f"-dump 0 -o {gro}'"
            )
        else:
            exec_str = (
                f"executing 'echo System | gmx trjconv -s {top} -f {traj}"
                f" -dump 0 -o {gro}'"
            )
        logger.debug(exec_str)
        os.system(exec_str)
        try:
            u = Universe(gro, traj)
            # write first frame into gro file
            u.atoms.write(gro, frames=u.trajectory[[0]])
        except Exception as e:
            logger.warning(e)
        struc = gro

    # if there is a topology and MDAnalysis reads it, we can use zero-frame
    # extraction as a structure!
    if not fail_from_top and struc is None:
        struc = gro

    # At last, we create universe from just a structure!
    logger.info(f"Making Universe from {struc}..")
    u0 = Universe(struc)

    lipids = []

    # select lipids
    for key_mol in lipids_set.names:
        logger.info(f"Calculating number of '{key_mol}' lipids")
        selection = ""
        if key_mol in sim["COMPOSITION"]:
            lip = Lipid(key_mol)
            m_file = sim["COMPOSITION"][key_mol]["MAPPING"]
            lip.register_mapping(m_file)
            for key in lip.mapping_dict:
                if "RESIDUE" in lip.mapping_dict[key]:
                    selection = (
                        selection
                        + "(resname "
                        + lip.mapping_dict[key]["RESIDUE"]
                        + " and name "
                        + lip.mapping_dict[key]["ATOMNAME"]
                        + ") or "
                    )
                else:
                    selection = "resname " + sim["COMPOSITION"][key_mol]["NAME"]
                    break
        selection = selection.rstrip(" or ")
        molecules = u0.select_atoms(selection)
        if molecules.n_residues > 0:
            lipids.append(u0.select_atoms(selection))

    # join all the selected the lipids together to make a selection of the entire
    # membrane and calculate the z component of the centre of mass of
    # the membrane
    membrane = u0.select_atoms("")
    R_membrane_z = 0
    if not lipids:
        for i in range(0, len(lipids)):
            membrane = membrane + lipids[i]
        R_membrane_z = membrane.center_of_mass()[2]
    logger.info(f"Center of the mass of the membrane: {str(R_membrane_z)}")

    # ---- number of each lipid per leaflet
    # TODO: remove code duplication block!
    for key_mol in lipids_set.names:
        leaflet1 = 0
        leaflet2 = 0

        selection = ""
        if key_mol in sim["COMPOSITION"]:
            lip = Lipid(key_mol)
            m_file = sim["COMPOSITION"][key_mol]["MAPPING"]
            lip.register_mapping(m_file)
            for key in lip.mapping_dict:
                if "RESIDUE" in lip.mapping_dict[key]:
                    selection = (
                        selection
                        + "resname "
                        + lip.mapping_dict[key]["RESIDUE"]
                        + " and name "
                        + lip.mapping_dict[key]["ATOMNAME"]
                        + " or "
                    )
                    break
                else:
                    selection = "resname " + sim["COMPOSITION"][key_mol]["NAME"]
                    break

        # if lipid was found then selection is not empty
        if selection != "":
            selection = selection.rstrip(" or ")
            logger.debug(f"Selection: `{selection}`")
            molecules = u0.select_atoms(selection)
            logger.debug("Resnames: " +
                         ", ".join(molecules.residues.resnames) +
                         " | ResIDs: " +
                         ", ".join(map(str, molecules.residues.resids))
                         )

            if molecules.n_residues > 0:
                for mol in molecules.residues:
                    R = mol.atoms.center_of_mass()

                    if R[2] - R_membrane_z > 0:
                        leaflet1 = leaflet1 + 1
                    elif R[2] - R_membrane_z < 0:
                        leaflet2 = leaflet2 + 1

        try:
            sim["COMPOSITION"][key_mol]["COUNT"] = [leaflet1, leaflet2]
        except KeyError:
            continue
        else:
            logger.info(f"Number of '{key_mol}' in upper leaflet: {str(leaflet1)}")
            logger.info(f"Number of '{key_mol}' in lower leaflet: {str(leaflet2)}")

    # ----- numbers of other molecules

    for key_mol in molecules_set.names:
        try:
            mol_name = sim["COMPOSITION"][key_mol]["NAME"]
        except KeyError:
            continue
        else:
            mol_number = u0.select_atoms("resname " + mol_name).n_residues
            sim["COMPOSITION"][key_mol]["COUNT"] = mol_number
            logger.info(
                f"Number of '{key_mol}': {str(sim['COMPOSITION'][key_mol]['COUNT'])}"
            )

    # Anne: Read trajectory size and length

    sim["TRAJECTORY_SIZE"] = os.path.getsize(traj)

    dt = 0
    nsteps = 0
    nstxout = 0

    n_frames = len(u.trajectory)
    timestep = u.trajectory.dt

    logger.info(f"Number of frames: {n_frames}")
    logger.info(f"Timestep: {timestep}")

    trj_length = n_frames * timestep

    sim["TRJLENGTH"] = trj_length

    # Read temperature from tpr
    if sim["SOFTWARE"].upper() == "GROMACS":
        file1 = os.path.join(dir_tmp, "tpr.txt")

        logger.info("Exporting information with gmx dump")
        if (
            "WARNINGS" in sim and
            sim["WARNINGS"] is not None and
            "GROMACS_VERSION" in sim["WARNINGS"] and
            sim["WARNINGS"]["GROMACS_VERSION"] == "gromacs3"
           ):
            os.system(f"echo System | gmxdump -s {top} > {file1}")
            TemperatureKey = "ref_t"
        else:
            os.system(f"echo System | gmx dump -s {top} > {file1}")
            TemperatureKey = "ref-t"

        with open(file1, "rt") as tpr_info:
            for line in tpr_info:
                if TemperatureKey in line:
                    sim["TEMPERATURE"] = float(line.split()[1])

    logger.info("Parameters read from input files:")
    logger.info(f"TEMPERATURE: {str(sim['TEMPERATURE'])}")
    logger.info(f"LENGTH OF THE TRAJECTORY: {str(sim['TRJLENGTH'])}")

    # Check that the number of atoms between data and README.yaml match

    natoms_trj = u.atoms.n_atoms

    number_of_atoms = 0
    for key_mol in sim["COMPOSITION"]:
        mol = Lipid(key_mol) if key_mol in lipids_set else NonLipid(key_mol)
        mol.register_mapping(sim["COMPOSITION"][key_mol]["MAPPING"])

        if sim.get("UNITEDATOM_DICT") and "SOL" not in key_mol:
            mapping_file_length = 0

            for key in mol.mapping_dict:
                if "H" in key:
                    continue
                else:
                    mapping_file_length += 1

        else:
            mapping_file_length = len(mol.mapping_dict)

        number_of_atoms += (
                np.sum(sim["COMPOSITION"][key_mol]["COUNT"]) * mapping_file_length
            )

    if number_of_atoms != natoms_trj:
        stop = input(
            f"Number of atoms in trajectory {natoms_trj} and README.yaml "
            f"{number_of_atoms} do no match. Check the mapping files and molecule"
            f" names. {os.linesep} If you know what you are doing, you can still "
            "continue the running the script. Do you want to (y/n)?"
        )
        if stop == "n":
            os._exit("Interrupted because atomnumbers did not match")
        if stop == "y":
            logger.warning(
                "Progressed even thought that atom numbers did not match."
                " CHECK RESULTS MANUALLY!"
            )

    sim["NUMBER_OF_ATOMS"] = natoms_trj
    logger.info(f"Number of atoms in the system: {str(sim['NUMBER_OF_ATOMS'])}")

    # ---- DATE OF RUNNING ----
    today = date.today().strftime("%d/%m/%Y")
    sim["DATEOFRUNNING"] = today

    logger.info(f"Date of adding to the databank: {sim['DATEOFRUNNING']}")

    # Type of system is currently hard coded because only lipid bilayers are currently
    # added. When we go for other systems, this will be given by user.
    if "TYPEOFSYSTEM" not in list(sim.keys()):
        sim["TYPEOFSYSTEM"] = "lipid bilayer"

    # ---- Save to databank

    try:
        directory_path = create_databank_directories(sim, sim_hashes, args.output_dir)
    except NotImplementedError as e:
        logger.error(e)
        quit(4)
    except OSError as e:
        logger.error(f"couldn't create output directory: {e.args[1]}")
        quit(2)

    logger.info(f"saving results to '{directory_path}'")

    # copy previously downloaded files
    logger.info("copying previously downloaded files ...")
    shutil.copyfile(traj, os.path.join(directory_path, os.path.basename(traj)))
    shutil.copyfile(top, os.path.join(directory_path, os.path.basename(top)))

    # dictionary saved in yaml format
    outfile_dict = os.path.join(dir_tmp, "README.yaml")

    with open(outfile_dict, "w") as f:
        yaml.dump(sim, f, sort_keys=False, allow_unicode=True)
        shutil.copyfile(
            os.path.join(dir_tmp, "README.yaml"),
            os.path.join(directory_path, "README.yaml"),
        )

    logger.info("Script completed successfully!")
