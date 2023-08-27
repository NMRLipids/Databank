# Library that contains all dictionaries and functions used in building and analyzing the NMR lipids databank

####################DICTIONARIES#####################################
# Define dictionaries
#
# Dictionary of lipids.
#
# If you add a lipid which is not yet in the databank, you have to add it here

import copy
import hashlib
from io import BufferedReader
import socket, urllib, logging
import shutil
from urllib.error import HTTPError

import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path
from typing import List

logger = logging.getLogger("__name__")

lipids_dict = {
    "POPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "POPG": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "POPS": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "POPE": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "PYPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "PAzePCprot": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "PAzePCdeprot": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DMPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DPPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DPPE": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DPPG": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DEPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DRPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DYPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DLPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DLIPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DOG": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DOPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DOPE": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DDOPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DOPS": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DSPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DAPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DMTAP": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SDG": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SDPE": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SOPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "POPI": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SAPI": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SLPI": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "CER": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "CHOL": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DCHOL": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DHMDMAB": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DPPG": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SLiPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
}


# Dictionary of other than lipid molecules.
#
# If you add other than a lipid molecule which is not yet in the databank, you have to add it here

molecules_dict = {
    "POT": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SOD": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "CLA": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "CAL": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "CES": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SOL": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "C20": {
        "REQUIRED": False,
        "TYPE": "string",
    },
}


# Dictionary containing the force fields for molecules given by the contributor

molecule_ff_dict = {
    "FFPOPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFPOPG": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFPOPS": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFPOPE": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFPYPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFPAzePCprot": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFPAzePCdeprot": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDPPG": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDMPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDPPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDPPE": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDEPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDRPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDYPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDLPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDLIPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDOG": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDOPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDOPE": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDDOPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDOPS": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDSPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDAPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDMTAP": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFSOPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFPOPI": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFSDG": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFSDPE": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFSAPI": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFSLPI": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFCER": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFCHOL": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDCHOL": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDHMDMAB": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDPPG": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFPOT": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFSOD": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFCLA": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFCAL": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFCES": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFSOL": {
        "REQUIRED": False,
        "TYPE": "string",
    },
}


# Databank dictionary for simulations ran with Gromacs

gromacs_dict = {
    "INI": {
        "REQUIRED": False,
        "TYPE": "files",
        "EXTENSION": (
            "gro",
            "pdb",
        ),
    },  # Could be not needed in the future (tpr)
    "MDP": {
        "REQUIRED": False,
        "TYPE": "file",
        "EXTENSION": ("mdp",),
    },  # Could be not needed in the future (tpr)
    "TRJ": {
        "REQUIRED": True,
        "TYPE": "files",
        "EXTENSION": (
            "xtc",
            "trr",
        ),
    },
    "TPR": {
        "REQUIRED": True,
        "TYPE": "file",
        "EXTENSION": ("tpr",),
    },
    "CPT": {
        "REQUIRED": False,
        "TYPE": "file",
        "EXTENSION": ("cpt",),
    },
    "TOP": {
        "REQUIRED": False,
        "TYPE": "file",
        "EXTENSION": ("top",),
    },
    "ITP": {
        "REQUIRED": False,
        "TYPE": "files",
        "EXTENSION": ("itp",),
    },
    "LOG": {
        "REQUIRED": False,
        "TYPE": "file",
        "EXTENSION": ("log",),
    },
    "GRO": {
        "REQUIRED": False,
        "TYPE": "file",
        "EXTENSION": ("gro",),
    },
    "FF": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FF_SOURCE": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FF_DATE": {
        "REQUIRED": False,
        "TYPE": "date",
    },
    "DOI": {
        "REQUIRED": True,
        "TYPE": "string",
    },
    "SYSTEM": {
        "REQUIRED": True,
        "TYPE": "string",
    },
    "TEMPERATURE": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "TRJLENGTH": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "PREEQTIME": {
        "REQUIRED": True,
        "TYPE": "integer",
    },
    "TIMELEFTOUT": {
        "REQUIRED": True,
        "TYPE": "integer",
    },
    "UNITEDATOM_DICT": {
        "REQUIRED": False,
        "TYPE": "dictionary",
    },
    "PUBLICATION": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "AUTHORS_CONTACT": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SOFTWARE_VERSION": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DATEOFRUNNING": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "NUMBER_OF_ATOMS": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "TRAJECTORY_SIZE": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "DIR_WRK": {
        "REQUIRED": True,
        "TYPE": "string",
    },
    "COMPOSITION": {
        "REQUIRED": True,
        "TYPE": "dictionary",
    },
    "WARNINGS": {
        "REQUIRED": False,
        "TYPE": "dictionary",
    },
}

# Amber
amber_dict = {
    "TRJ": {
        "REQUIRED": True,
        "TYPE": "files",
        "EXTENSION": (
            "nc",
            "ncdf",
            "trj",
            "mdcrd",
        ),
    },
    "TOP": {
        "REQUIRED": False,
        "TYPE": "file",
        "EXTENSION": (
            "prmtop",
            "top",
            "parm7",
        ),
    },
    "FF": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FF_SOURCE": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FF_DATE": {
        "REQUIRED": False,
        "TYPE": "date",
    },
    "PDB": {
        "REQUIRED": True,
        "TYPE": "file",
        "EXTENSION": "pdb",
    },
    "DOI": {
        "REQUIRED": True,
        "TYPE": "string",
    },
    "SYSTEM": {
        "REQUIRED": True,
        "TYPE": "string",
    },
    "TEMPERATURE": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "TRJLENGTH": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "PREEQTIME": {
        "REQUIRED": True,
        "TYPE": "integer",
    },
    "TIMELEFTOUT": {
        "REQUIRED": True,
        "TYPE": "integer",
    },
    "UNITEDATOM_DICT": {
        "REQUIRED": False,
        "TYPE": "dictionary",
    },
    "PUBLICATION": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "AUTHORS_CONTACT": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SOFTWARE_VERSION": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DATEOFRUNNING": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "NUMBER_OF_ATOMS": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "TRAJECTORY_SIZE": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "DIR_WRK": {
        "REQUIRED": True,
        "TYPE": "string",
    },
    "COMPOSITION": {
        "REQUIRED": True,
        "TYPE": "dictionary",
    },
}

# NAMD
namd_dict = {
    "TRJ": {
        "REQUIRED": True,
        "TYPE": "files",
        "EXTENSION": ("dcd"),
    },
    "INP": {
        "REQUIRED": False,
        "TYPE": "file",
        "EXTENSION": (".inp"),
    },
    "LOG": {
        "REQUIRED": False,
        "TYPE": "files",
        "EXTENSION": ("log"),
        # can be parsed to get software version etc.
    },
    "TOP": {
        "REQUIRED": False,
        "TYPE": "file",
        "EXTENSION": ("psf"),
    },
    "FF": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FF_SOURCE": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FF_DATE": {
        "REQUIRED": False,
        "TYPE": "date",
    },
    "PDB": {
        "REQUIRED": True,
        "TYPE": "file",
        "EXTENSION": "pdb",
    },
    "DOI": {
        "REQUIRED": True,
        "TYPE": "string",
    },
    "SYSTEM": {
        "REQUIRED": True,
        "TYPE": "string",
    },
    "TEMPERATURE": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "TRJLENGTH": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "PREEQTIME": {
        "REQUIRED": True,
        "TYPE": "integer",
    },
    "TIMELEFTOUT": {
        "REQUIRED": True,
        "TYPE": "integer",
    },
    "UNITEDATOM_DICT": {
        "REQUIRED": False,
        "TYPE": "dictionary",
    },
    "PUBLICATION": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "AUTHORS_CONTACT": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SOFTWARE_VERSION": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DATEOFRUNNING": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "NUMBER_OF_ATOMS": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "TRAJECTORY_SIZE": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "DIR_WRK": {
        "REQUIRED": True,
        "TYPE": "string",
    },
    "COMPOSITION": {
        "REQUIRED": True,
        "TYPE": "dictionary",
    },
}

# CHARMM
charmm_dict = {}

# OPENMM
openmm_dict = {
    "TRJ": {
        "REQUIRED": True,
        "TYPE": "files",
        "EXTENSION": ("xtc", "trr", "nc", "ncdf", "trj", "mdcrd", "dcd"),
    },
    "PDB": {
        "REQUIRED": True,
        "TYPE": "file",
        "EXTENSION": ("pdb",),
    },
    "TOP": {
        "REQUIRED": False,
        "TYPE": "file",
        "EXTENSION": ("psf",),
    },
    "XML": {
        "REQUIRED": False,  # state files from openmm, almost similar to a restart file
        "TYPE": "file",
        "EXTENSION": ("xml",),
    },
    "CHK": {
        "REQUIRED": False,
        "TYPE": "file",
        "EXTENSION": ("chk",),
    },
    "CRD": {
        "REQUIRED": False,
        "TYPE": "file",
        "EXTENSION": ("crd",),
    },
    "INP": {
        "REQUIRED": False,  # input file used to run the simulation
        "TYPE": "file",
        "EXTENSION": ("inp",),
    },
    "FF": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FF_SOURCE": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FF_DATE": {
        "REQUIRED": False,
        "TYPE": "date",
    },
    "DOI": {
        "REQUIRED": True,
        "TYPE": "string",
    },
    "SYSTEM": {
        "REQUIRED": True,
        "TYPE": "string",
    },
    "TEMPERATURE": {
        "REQUIRED": False,
        "TYPE": "float",
    },
    "TRJLENGTH": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "PREEQTIME": {
        "REQUIRED": True,
        "TYPE": "integer",
    },
    "TIMELEFTOUT": {
        "REQUIRED": True,
        "TYPE": "integer",
    },
    "UNITEDATOM_DICT": {
        "REQUIRED": False,
        "TYPE": "dictionary",
    },
    "PUBLICATION": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "AUTHORS_CONTACT": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SOFTWARE_VERSION": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DATEOFRUNNING": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "NUMBER_OF_ATOMS": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "TRAJECTORY_SIZE": {
        "REQUIRED": False,
        "TYPE": "integer",
    },
    "DIR_WRK": {
        "REQUIRED": True,
        "TYPE": "string",
    },
    "COMPOSITION": {
        "REQUIRED": True,
        "TYPE": "dictionary",
    },
}

# SOFTWARE
software_dict = {
    "GROMACS": gromacs_dict,
    "AMBER": amber_dict,
    "NAMD": namd_dict,
    "CHARMM": charmm_dict,
    "OPENMM": openmm_dict,
}


##############CLASS FOR LOOPING OVER SYSTEMS#######################################

import yaml


class databank:
    """ After intialized, this class contains the list of simulations in the NMR lipids databank. 
    
    As an argument, it takes the local location of /Data/Simulations/ folder.
    
    Example usage to loop over systems: 

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

    def __load_systems__(self, path):
        for subdir, dirs, files in os.walk(path):
            for filename in files:
                filepath = os.path.join(subdir, filename)
                # print(filepath)
                if filename == "README.yaml":
                    with open(filepath) as yaml_file:
                        content = yaml.load(yaml_file, Loader=yaml.FullLoader)
                        size = len(filepath)
                        content["path"] = filepath[: size - 11]
                        self.systems.append(content)

    def get_systems(self):
        return self.systems

    def pie_temperature(self):
        list_feature = [int(float(system["TEMPERATURE"])) for system in self.systems]
        import collections

        counter = collections.Counter(list_feature)
        plt.pie(counter.values(), labels=counter.keys(), normalize=True)


#########################FUNCTIONS###################################################
# functions used in building and analyzing the databank


# Download link
def download_link(doi, file):  # deprecated?
    if "zenodo" in doi.lower():
        zenodo_entry_number = doi.split(".")[2]
        return "https://zenodo.org/record/" + zenodo_entry_number + "/files/" + file
    else:
        print("DOI provided: {0}".format(doi))
        print(
            "Repository not validated. Please upload the data for example to zenodo.org"
        )
        return ""


def resolve_doi_url(doi: str, validate_uri: bool = True) -> str:
    """Returns full doi link of given ressource, also checks if URL is valid.

    Args:
        doi (str): [doi] part from config
        validate_uri (bool, optional): Check if URL is valid. Defaults to True.

    Returns:
        str: full doi link
    """
    res = "https://doi.org/" + doi

    if validate_uri:
        socket.setdefaulttimeout(10)  # seconds
        _ = urllib.request.urlopen(res)
    return res


def resolve_download_file_url(doi: str, fi_name: str, validate_uri: bool = True) -> str:
    """Resolve file URI from supported DOI with given filename

    Args:
        doi (str): DOI string
        fi_name (str): name of the file to resolve from source
        validate_uri (bool, optional): Check if URI exists. Defaults to True.

    Raises:
        NotImplementedError: Unsupported DOI repository
        HTTPError: HTTP Error Status Code
        URLError: Failed to reach the server

    Returns:
        str: file URI
    """
    if "zenodo" in doi.lower():
        zenodo_entry_number = doi.split(".")[2]
        uri = "https://zenodo.org/record/" + zenodo_entry_number + "/files/" + fi_name

        # check if ressource exists, may throw exception
        if validate_uri:
            socket.setdefaulttimeout(10)  # seconds
            res = urllib.request.urlopen(uri)
        return uri
    else:
        raise NotImplementedError(
            "Repository not validated. Please upload the data for example to zenodo.org"
        )


def download_resource_from_uri(
    uri: str, dest: str, override_if_exists: bool = False
) -> None:
    """Download file resource [from uri] to given file destination using urllib

    Args:
        uri (str): file URL
        dest (str): file destination path
        override_if_exists (bool, optional): Override dest. file if exists. Defaults to False.

    Raises:
        Exception: HTTPException: An error occured during download

    Returns:
        None
    """
    # TODO verify file size before skipping already existing download!

    class RetrieveProgressBar(tqdm):
        # uses tqdm.update(), see docs https://github.com/tqdm/tqdm#hooks-and-callbacks
        def update_retrieve(self, b=1, bsize=1, tsize=None):
            if tsize is not None:
                self.total = tsize
            return self.update(b * bsize - self.n)

    fi_name = uri.split("/")[-1]

    # check if dest path already exists
    if not override_if_exists and os.path.isfile(dest):
        socket.setdefaulttimeout(10)  # seconds

        # compare filesize
        fi_size = urllib.request.urlopen(uri).length  # download size
        if fi_size == os.path.getsize(dest):
            logger.info(f"{dest}: file already exists, skipping")
            return
        else:
            logger.warning(
                f"{fi_name} filesize mismatch of local file '{fi_name}', redownloading ..."
            )

    # download
    socket.setdefaulttimeout(10)  # seconds

    url_size = urllib.request.urlopen(uri).length  # download size

    with RetrieveProgressBar(
        unit="B", unit_scale=True, unit_divisor=1024, miniters=1, desc=fi_name
    ) as u:
        response = urllib.request.urlretrieve(uri, dest, reporthook=u.update_retrieve)

    # check if the file is fully downloaded
    size = os.path.getsize(dest)

    if url_size != size:
        raise Exception(f"downloaded filsize mismatch ({size}/{url_size} B)")


class YamlBadConfigException(Exception):
    """Custom Exception class for parsing the yaml configuration"""

    def __init__(self, *args, **kwargs) -> None:
        Exception.__init__(self, *args, **kwargs)


def parse_valid_config_settings(info_yaml: dict) -> (dict, List[str]):
    """Parses, validates and updates dict entries from yaml configuration file.

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
        missing_keys = [k for k in list(sim.keys()) if k not in software_required_keys]
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


# Return mapping name of atom from mapping file
def read_mapping_file(mapping_file, atom1):
    with open(mapping_file, "rt") as mapping_file:
        mapping = yaml.load(mapping_file, Loader=yaml.FullLoader)
        m_atom1 = mapping[atom1]["ATOMNAME"]
    return m_atom1


# Return mapping names of pair of atoms from mapping file
def read_mapping_filePAIR(mapping_file, atom1, atom2):
    with open(mapping_file, "rt") as mapping_file:
        print(mapping_file)
        for line in mapping_file:
            if atom1 in line:
                m_atom1 = line.split()[1]
            #                    print(m_atom1)
            if atom2 in line:
                m_atom2 = line.split()[1]
    #                    print(m_atom2)
    return m_atom1, m_atom2


def read_mapping_file_res(mapping_file, atom1):
    with open(mapping_file, "rt") as mapping_file:
        for line in mapping_file:
            if atom1 in line:
                m_res = line.split()[2]
    return m_res


def make_positive_angles(x):
    for i in range(len(x)):
        if x[i] < 0:
            x[i] = x[i] + 360
        else:
            x[i] = x[i]
    return x


######################DIHEDRALS##############################################
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import calc_dihedrals
from MDAnalysis.analysis.dihedrals import Dihedral


class DihedralFromAtoms(AnalysisBase):
    """Calculate dihedral angles for specified atomgroups.

    Dihedral angles will be calculated for each atomgroup that is given for
    each step in the trajectory. Each :class:`~MDAnalysis.core.groups.AtomGroup`
    must contain 4 atoms.

    Note
    ----
    This class takes a list as an input and is most useful for a large
    selection of atomgroups. If there is only one atomgroup of interest, then
    it must be given as a list of one atomgroup.

    """

    def __init__(self, atomgroups, orders, **kwargs):
        """Parameters
        ----------
        atomgroups : list
                        a list of atomgroups for which the dihedral angles are calculated

        Raises
        ------
        ValueError
                        If any atomgroups do not contain 4 atoms

        """
        super(DihedralFromAtoms, self).__init__(
            atomgroups[0].universe.trajectory, **kwargs
        )
        self.atomgroups = atomgroups

        if any([len(ag) != 4 for ag in atomgroups]):
            raise ValueError("All AtomGroups must contain 4 atoms")

        if len(atomgroups) != len(orders):
            raise ValueError("Order data should be provided for every atom group")

        self.ag1 = mda.AtomGroup(
            [atomgroups[i][orders[i][0]] for i in range(len(atomgroups))]
        )
        self.ag2 = mda.AtomGroup(
            [atomgroups[i][orders[i][1]] for i in range(len(atomgroups))]
        )
        self.ag3 = mda.AtomGroup(
            [atomgroups[i][orders[i][2]] for i in range(len(atomgroups))]
        )
        self.ag4 = mda.AtomGroup(
            [atomgroups[i][orders[i][3]] for i in range(len(atomgroups))]
        )

    def _prepare(self):
        self.angles = []

    def _single_frame(self):
        angle = calc_dihedrals(
            self.ag1.positions,
            self.ag2.positions,
            self.ag3.positions,
            self.ag4.positions,
            box=self.ag1.dimensions,
        )
        self.angles.append(angle)

    def _conclude(self):
        self.angles = np.rad2deg(np.array(self.angles))


def calcDihedrals(lipids, DIHatoms):
    colors = {"POPC": "black", "POPS": "red", "POPE": "blue", "POPG": "green"}

    found = 0

    for subdir, dirs, files in os.walk(r"../../Data/Simulations/"):
        for filename in files:
            if filename.endswith("README.yaml"):
                READMEfilepath = subdir + "/README.yaml"
                found += 1
                print(READMEfilepath)
                with open(READMEfilepath) as yaml_file:
                    readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
                    for molname in lipids:
                        doi = readme.get("DOI")
                        trj = readme.get("TRJ")
                        tpr = readme.get("TPR")
                        trj_name = subdir + "/" + readme.get("TRJ")[0][0]
                        tpr_name = subdir + "/" + readme.get("TPR")[0][0]
                        gro_name = subdir + "/conf.gro"
                        trj_url = download_link(doi, trj[0][0])
                        tpr_url = download_link(doi, tpr[0][0])

                        # Download tpr and xtc files to same directory where dictionary and data are located
                        if not os.path.isfile(tpr_name):
                            response = urllib.request.urlretrieve(tpr_url, tpr_name)

                        if not os.path.isfile(trj_name):
                            response = urllib.request.urlretrieve(trj_url, trj_name)

                        if sum(readme["N" + molname]) > 0:
                            print("Analyzing " + molname + " in " + READMEfilepath)
                            # fig= plt.figure(figsize=(12,9))
                            if not os.path.isfile(gro_name):
                                os.system(
                                    "echo System | gmx trjconv -f {} -s {}  -dump 0 -o {}".format(
                                        trj_name, tpr_name, gro_name
                                    )
                                )

                            xtcwhole = subdir + "/whole.xtc"
                            if not os.path.isfile(xtcwhole):
                                os.system(
                                    "echo System | gmx trjconv -f {} -s {} -o {} -pbc mol ".format(
                                        trj_name, tpr_name, xtcwhole
                                    )
                                )

                            try:
                                traj = mda.Universe(tpr_name, xtcwhole)
                            except FileNotFoundError or OSError:
                                continue

                            mapping_file = (
                                "../BuildDatabank/mapping_files/"
                                + readme["MAPPING_DICT"][molname]
                            )  # readme.get('MAPPING')[0][0]
                            try:
                                atom1 = read_mapping_file(mapping_file, DIHatoms[0])
                                atom2 = read_mapping_file(mapping_file, DIHatoms[1])
                                atom3 = read_mapping_file(mapping_file, DIHatoms[2])
                                atom4 = read_mapping_file(mapping_file, DIHatoms[3])
                                print(atom1, atom2, atom3, atom4)
                            except:
                                print("Some atom not found in the mapping file.")
                                continue

                            ags = []
                            orders = []
                            for residue in traj.select_atoms(
                                "name {} {} {} {}".format(atom1, atom2, atom3, atom4)
                            ).residues:
                                atoms = traj.select_atoms(
                                    "name {} {} {} {} and resid {}".format(
                                        atom1, atom2, atom3, atom4, str(residue.resid)
                                    )
                                )
                                if len(atoms) == 4:
                                    # print(atoms.names)
                                    ags.append([])
                                    ags[-1] = copy.deepcopy(atoms)
                                    orders.append(
                                        [
                                            np.where(atoms.names == atom1)[0][0],
                                            np.where(atoms.names == atom2)[0][0],
                                            np.where(atoms.names == atom3)[0][0],
                                            np.where(atoms.names == atom4)[0][0],
                                        ]
                                    )
                            R = DihedralFromAtoms(ags, orders).run()
                            dihRESULT = R.angles.T

                            dihRESULT = [make_positive_angles(x) for x in dihRESULT]
                            distSUM = np.zeros(360)
                            for i in dihRESULT:
                                distSUM += np.histogram(
                                    i, np.arange(0, 361, 1), density=True
                                )[0]

                            distSUM = [x / len(dihRESULT) for x in distSUM]
                            xaxis = (
                                np.histogram(i, np.arange(0, 361, 1), density=True)[1][
                                    1:
                                ]
                                + np.histogram(i, np.arange(0, 361, 1), density=True)[
                                    1
                                ][:-1]
                            ) / 2.0

                            dihedralFOLDERS = subdir.replace("Simulations", "dihedral")
                            os.system("mkdir -p {}".format(dihedralFOLDERS))
                            os.system(
                                "cp {} {}".format(READMEfilepath, dihedralFOLDERS)
                            )
                            outfile = open(
                                str(dihedralFOLDERS)
                                + "/"
                                + molname
                                + "_"
                                + DIHatoms[0]
                                + "_"
                                + DIHatoms[1]
                                + "_"
                                + DIHatoms[2]
                                + "_"
                                + DIHatoms[3]
                                + ".dat",
                                "w",
                            )
                            for i in range(len(xaxis)):
                                outfile.write(
                                    str(xaxis[i]) + " " + str(distSUM[i]) + "\n"
                                )
                            outfile.close()


#######################ORDER PARAMETERS######################################
"""
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
        # print(op.resname + " " + op.atAname + " " + op.atBname)
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
        Nres = len(op.selection)
        op.traj = [0] * Nres

    for frame in mol.trajectory:
        for op in ordPars:  # .values():
            Nres = len(op.selection)
            for i in range(0, Nres):
                residue = op.selection[i]
                S = op.calc_OP(residue)
                op.traj[i] = op.traj[i] + S / Nframes


def parse_op_input(mapping_file, lipid_name):
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
    regexp2_C = re.compile(r"M_G[0-9]_M")
    regexp3_C = re.compile(r"M_C[0-9]_M")

    for mapping_key in mapping_dict.keys():
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
            #            print(atomC)
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

        if atomH:
            items = [atomC[1], atomH[1], atomC[0], atomH[0]]
            #               print(resname + " " + atomH[1])
            op = OrderParameter(resname, items[0], items[1], items[2], items[3])
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
    ordPars = parse_op_input(inp_fname, lipid_name)

    read_trajs_calc_OPs(ordPars, top_fname, traj_fname)

    return ordPars


#######################################################################################


def calc_angle(atoms, com):
    """
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
    u = mda.Universe(gro)
    z = u.dimensions[2]
    return z


def read_trj_PN_angles(molname, atoms, top_fname, traj_fname, gro_fname):
    mol = mda.Universe(top_fname, traj_fname)

    selection = mol.select_atoms(
        "resname " + molname + " and (name " + atoms[0] + ")",
        "resname " + molname + " and (name " + atoms[1] + ")",
    ).atoms.split("residue")
    com = mol.select_atoms(
        "resname " + molname + " and (name " + atoms[0] + " or name " + atoms[1] + ")"
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
    """Calculates sha1 hash of given file using hashlib

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
    """create nested output directory structure to save results

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
    elif sim["SOFTWARE"] == "openMM":
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
