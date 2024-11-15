"""
Library contains all API functions and many functions used in building and
analyzing the NMRlipids databank
"""

import copy
import hashlib
import urllib
import logging
from tqdm import tqdm

import json
import yaml
import os
import sys
import numpy as np
import math
import MDAnalysis as mda

from DatabankLib import NMLDB_SIMU_PATH, NMLDB_ROOT_PATH
from DatabankLib.databank_defs import (
    lipids_dict, software_dict, molecules_dict, molecule_ff_dict)
from DatabankLib.databankio import resolve_download_file_url

logger = logging.getLogger(__name__)


def CalcAreaPerMolecule(system):
    """
    Calculates average area per lipid for a simulation defined with ``system``.
    It is using the ``apl.json`` file where area per lipid as a function of time
    calculated by the ``calcAPL.py`` is stored.

    :param system: NMRlipids databank dictionary defining a simulation.

    :return: area per lipid (Å^2)
    """
    APLpath = os.path.join(NMLDB_SIMU_PATH, system['path'], 'apl.json')
    try:
        with open(APLpath) as f:
            APLdata = json.load(f)
        sumAPL = 0
        sumIND = 0
        for i, j in APLdata.items():
            sumAPL += j
            sumIND += 1
        APL = sumAPL/sumIND
        return APL
    except Exception:
        print('apl.json not found from' + APLpath)


def GetThickness(system):
    """
    Gets thickness for a simulation defined with ``system`` from the ``thickness.json``
    file where thickness calculated by the ``calc_thickness.py`` is stored.

    :param system: NMRlipids databank dictionary defining a simulation.

    :return: membrane thickess (nm) or None
    """
    ThicknessPath = os.path.join(NMLDB_SIMU_PATH, system['path'], 'thickness.json')
    try:
        with open(ThicknessPath) as f:
            thickness = json.load(f)
        return thickness
    except Exception:
        return None


def ShowEquilibrationTimes(system: dict):
    """
    Prints relative equilibration time for each lipid within a simulation defined
    by ``system``. Relative equilibration times are calculated with
    ``NMRPCA_timerelax.py`` and stored in ``eq_times.json`` files.

    :param system: NMRlipids databank dictionary defining a simulation.
    """

    EqTimesPath = os.path.join(NMLDB_SIMU_PATH, system['path'], 'eq_times.json')

    try:
        with open(EqTimesPath) as f:
            EqTimeDict = json.load(f)
    except Exception:
        raise FileNotFoundError(f'eq_times.json not found for {system["ID"]}')

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
    Returns a string using MDAnalysis notation that can used to select all lipids from
    the ``system``.

    :param system: NMRlipids databank dictionary defining a simulation.

    :return: a string using MDAnalysis notation that can used to select all lipids from
             the ``system``.
    """

    resSet = set()
    for key in system['COMPOSITION'].keys():
        if key in molecules:
            m_file = system['COMPOSITION'][key]['MAPPING']
            mapping_dict = loadMappingFile(m_file)
            try:
                for atom in mapping_dict:
                    resSet.add(mapping_dict[atom]['RESIDUE'])
            except (KeyError, TypeError):
                resSet.add(system['COMPOSITION'][key]['NAME'])

    lipids = 'resname ' + ' or resname '.join(sorted(list(resSet)))

    return lipids


def simulation2universal_atomnames(system, molecule, atom):
    """
    Get force field specific atom name corresponding to universal atom name from
    the ``system``.

    :param mapping_file: path for the mapping file
    :param atom1: universal atom name

    :return: force field specific atom name
    """
    try:
        mapping = loadMappingFile(system['COMPOSITION'][molecule]['MAPPING'])
    except Exception:
        sys.stderr.write('Mapping file was not found!\n')
        return None

    try:
        m_atom1 = mapping[atom]["ATOMNAME"]
    except (KeyError, TypeError):
        sys.stderr.write(
            f"{atom} was not found from {system['COMPOSITION'][molecule]['MAPPING']}!")
        return None

    return m_atom1


def loadMappingFile(mapping_file):
    """
    Load mapping file into a dictionary

    :param: name of the mapping file

    :return: mapping dictionary
    """
    mapping_file_path = os.path.join(
        NMLDB_ROOT_PATH, 'Scripts', 'DatabankLib', 'mapping_files', mapping_file)
    mapping_dict = {}
    with open(mapping_file_path, "r") as yaml_file:
        mapping_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)
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


def getUniversalAtomName(system: dict, atomName: str, lipid: str):
    """
    Returns the universal atom name corresponding the simulation specific ``atomName``
    of a ``lipid`` in a simulation defined by the ``system``.

    :param system: system dictionary
    :param atomName: simulation specific atomname
    :param lipid: universal lipid name

    :return: universal atomname (string) or None
    """
    try:
        mappingFile = system['COMPOSITION'][lipid]['MAPPING']
    except (KeyError, TypeError):
        sys.stderr.write('Mapping file was not found!\n')
        return None

    mappingDict = loadMappingFile(mappingFile)

    for universalName in mappingDict:
        simName = mappingDict[universalName]['ATOMNAME']
        if simName == atomName:
            return universalName

    sys.stderr.write('Atom was not found!\n')
    return None


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
    Takes the ``system`` dictionary as an input, downloads the required files to
    the NMRlipids databank directory and retuns MDAnalysis universe corressponding
    the ``system``.

    :param system: NMRlipids databank dictionary describing the simulation.

    :return: MDAnalysis universe
    """
    systemPath = os.path.join(NMLDB_SIMU_PATH, system['path'])
    doi = system.get('DOI')
    skipDownloading: bool = (doi == 'localhost')
    if skipDownloading:
        print("NOTE: The system with 'localhost' DOI should be downloaded by the user.")

    trj = system.get('TRJ')
    trj_name = os.path.join(systemPath, system.get('TRJ')[0][0])
    software = system['SOFTWARE']

    if (skipDownloading):
        if (not os.path.isfile(trj_name)):
            raise FileNotFoundError(
                f"Trajectory should be downloaded [{trj_name}] by user")
    else:
        trj_url = resolve_download_file_url(doi, trj[0][0])
        if (not os.path.isfile(trj_name)):
            print('Downloading trajectory with the size of ', system['TRAJECTORY_SIZE'],
                  ' to ', system['path'])
            _ = urllib.request.urlretrieve(trj_url, trj_name)

    if 'gromacs' in software:
        tpr_name = 'stub'
        try:
            tpr = system.get('TPR')
            tpr_name = os.path.join(systemPath, tpr[0][0])
            if skipDownloading:
                if (not os.path.isfile(tpr_name)):
                    raise FileNotFoundError(f"TPR should be downloaded [{tpr_name}]")
            else:
                tpr_url = resolve_download_file_url(doi, tpr[0][0])
                if (not os.path.isfile(tpr_name)):
                    _ = urllib.request.urlretrieve(tpr_url, tpr_name)

            u = mda.Universe(tpr_name, trj_name)
        except Exception:
            try:
                gro = system.get('GRO')
                conf = os.path.join(systemPath, gro[0][0])

                if skipDownloading:
                    if (not os.path.isfile(conf)):
                        raise FileNotFoundError(f"GRO should be downloaded [{conf}]")
                else:
                    gro_url = resolve_download_file_url(doi, gro[0][0])
                    if (not os.path.isfile(conf)):
                        _ = urllib.request.urlretrieve(gro_url, conf)
            except Exception:
                conf = os.path.join(systemPath, 'conf.gro')
            if (os.path.isfile(tpr_name)):
                print("Generating conf.gro because MDAnalysis cannot "
                      "(probably!) read tpr version")
                if (
                    'WARNINGS' in system and
                    'GROMACS_VERSION' in system['WARNINGS'] and
                    system['WARNINGS']['GROMACS_VERSION'] == 'gromacs3'
                ):
                    os.system(f'echo System | editconf -f {tpr_name} -o {conf}')
                else:
                    os.system(f"echo System | gmx trjconv "
                              f"-s {tpr_name} -f {trj_name} -dump 0 -o {conf}")
            u = mda.Universe(conf, trj_name)

    elif 'openMM' or 'NAMD' in software:
        pdb = system.get('PDB')
        pdb_name = os.path.join(systemPath, pdb[0][0])
        if skipDownloading:
            if (not os.path.isfile(pdb_name)):
                raise FileNotFoundError(f"PDB should be downloaded [{pdb_name}]")
        else:
            pdb_url = resolve_download_file_url(doi, pdb[0][0])
            if (not os.path.isfile(pdb_name)):
                _ = urllib.request.urlretrieve(pdb_url, pdb_name)
        u = mda.Universe(pdb_name, trj_name)

    else:
        raise NotImplementedError(
            'Other than GROMACS, openMM or NAMD are yet to be implemented.')

    return u


def read_trj_PN_angles(molname: str, atom1: str, atom2: str, MDAuniverse: mda.Universe):
    """
    Calculates the P-N vector angles with respect to membrane normal from the
    simulation defined by the MDAnalysis universe.

    :param molname: residue name of the molecule for which the P-N vector angle will
                    be calculated
    :param atom1: name of the P atom in the simulation
    :param atom2: name of the N atom in the simulation
    :param MDAuniverse: MDAnalysis universe of the simulation to be analyzed

    :return: tuple (angles of all molecules as a function of time,
                    time averages for each molecule,
                    the average angle over time and molecules,
                    the error of the mean calculated over molecules)
    """
    mol = MDAuniverse
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


# -------------------------------------- SEPARATED PART (??) ----------------------

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


def parse_valid_config_settings(info_yaml: dict) -> tuple[dict, list[str]]:
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
        missing_keys = [k for k in software_required_keys if k not in list(sim.keys())]
        raise YamlBadConfigException(
            f"Required '{sim['SOFTWARE'].upper()}' sim keys missing or "
            f"not defined in conf file: {', '.join(missing_keys)}"
        )

    logger.debug(
        f"all {len(software_required_keys)} required"
        f" '{sim['SOFTWARE'].upper()}' sim keys are present"
    )

    # STEP 3 - check working directory
    if "DIR_WRK" not in sim:
        raise KeyError("'DIR_WRK' Parameter missing in yaml")

    # STEP 4 - Check that all entry keys provided for each simulation are valid
    files_tbd = []

    #   loop config entries
    for key_sim, value_sim in sim.items():
        logger.debug(f"processing entry: sim['{key_sim}'] = {str(value_sim)}")

        if key_sim.upper() in "SOFTWARE":  # skip 'SOFTWARE' entry
            continue

        # STEP 4.1.
        # Anne: check if key is in molecules_dict, molecule_numbers_dict or
        # molecule_ff_dict too
        if (
            (key_sim.upper() not in software_sim.keys())
            and (key_sim.upper() not in molecules_dict.keys())
            and (key_sim.upper() not in lipids_dict.keys())
            and (key_sim.upper() not in molecule_ff_dict.keys())
        ):
            logger.error(
                f"key_sim '{key_sim}' in {sim['SOFTWARE'].lower()}_dict' "
                f": {key_sim.upper() in software_sim.keys()}"
            )
            logger.error(
                f"key_sim '{key_sim}' in molecules_dict "
                f": {key_sim.upper() in molecules_dict.keys()}"
            )
            logger.error(
                f"key_sim '{key_sim}' in lipids_dict "
                f": {key_sim.upper() in lipids_dict.keys()}"
            )
            logger.error(
                f"key_sim '{key_sim}' in molecule_ff_dict "
                f": {key_sim.upper() in molecule_ff_dict.keys()}"
            )
            raise YamlBadConfigException(
                f"'{key_sim}' not supported: Not found in "
                f"'{sim['SOFTWARE'].lower()}_dict', 'molecules_dict',"
                f" 'lipids_dict' and 'molecule_ff_dict'"
            )
        elif (
            key_sim.upper() not in software_sim.keys()
        ):  # hotfix for unkown yaml keys. TODO improve check 4.1?
            logger.warning(
                f"ignoring yaml entry '{key_sim}', not found "
                f"in '{sim['SOFTWARE'].lower()}_dict'"
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
                            f"found no file to download for "
                            f"entry '{key_sim}:{software_sim[key_sim]}'"
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
                # Batuhan: In conf file only one psf/tpr/pdb file allowed each
                # (can coexist), multiple TRJ files are ok
                # TODO true for all sim software?
                # TODO add dict entry param "unique" instead?
                if key_sim.upper() in ["PSF", "TPR", "PDB"] and len(sim[key_sim]) > 1:
                    raise YamlBadConfigException(
                        f"only one '{key_sim}' entry file allowed,"
                        f" but got {len(sim[key_sim])}: {sim[key_sim]}"
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
    formFactorPath = os.path.join(NMLDB_SIMU_PATH, system['path'], 'FormFactor.json')
    with open(formFactorPath) as f:
        formFactor = json.load(f)
    iprev = formFactor[0][1]
    iprevD = 0
    minX = []
    for i in formFactor:
        iD = i[1]-iprev
        if iD > 0 and iprevD < 0 and i[0] > 0.1:
            minX.append(i[0])
        iprevD = i[1] - iprev
        iprev = i[1]

    return minX


def averageOrderParameters(system):
    """
    Returns average order paramaters of *sn*-1 and *sn*-2 acyl chains based on universal
    atom names. The names starting with M_G1C will be assigned to sn-1 and names
    starting M_G2C to *sn*-2.

    :parameters system: a system dictionary

    :return: average of *sn*-1 and *sn*-2 order parameters
    """

    path = os.path.join(NMLDB_SIMU_PATH, system['path'])

    sn1sum = 0
    sn1count = 0
    sn2sum = 0
    sn2count = 0

    for lipid in system['COMPOSITION']:
        if lipid in lipids_dict and 'CHOL' not in lipid:
            OPpathSIM = os.path.join(path, lipid + 'OrderParameters.json')
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
    Returns hydration level of the system, i.e., number of water molecules divided
    by number of lipid molecules.

    :param system: a system dictionary

    :return: number of water molecules divided by number of lipid molecules
    """
    Nlipid = 0
    for molecule in system['COMPOSITION']:
        if molecule in lipids_dict:
            Nlipid += np.sum(system['COMPOSITION'][molecule]['COUNT'])
    Nwater = system['COMPOSITION']['SOL']['COUNT']
    return Nwater/Nlipid
