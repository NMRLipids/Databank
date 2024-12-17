"""
:module: settings/engines.py

:description: The module defines dictionaries describing parsing of YAML depending
              on MD engine used.

Engine-specific dictionaries have a fixed subfield types:
- REQUIRED
    Boolean field indicating if the field SHOULD be filled in `README.yaml`
- TYPE
    + file
    + files
    + string
    + dictionary
    + float
    + integer
- CATEGORY
    For file(s) type it can be:
    + structure
    + topology
    + trajectory
    + energy
    + None
- EXTENSION
    Is a tuple of possible file extensions
"""

import os.path
from collections.abc import Sequence
from typing import Optional

# GROMACS
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
        "CATEGORY": "trajectory",
        "EXTENSION": ("xtc", "trr",),
    },
    "TPR": {
        "REQUIRED": True,
        "TYPE": "file",
        "CATEGORY": "topology",
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
        "CATEGORY": "structure",
        "EXTENSION": ("gro",),
    },
    "EDR": {
        "REQUIRED": False,
        "TYPE": "file",
        "CATEGORY": "energy",
        "EXTENSION": ("edr",),
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
    "TYPEOFSYSTEM": {
        "REQUIRED": False,
        "TYPE": "string"
    },
    "BATCHID": {
        "REQUIRED": False,
        "TYPE": "string"
    }
}

# Amber
amber_dict = {
    "TRJ": {
        "REQUIRED": True,
        "TYPE": "files",
        "CATEGORY": "trajectory",
        "EXTENSION": ("nc", "ncdf", "trj", "mdcrd",),
    },
    "TOP": {
        "REQUIRED": False,
        "CATEGORY": "topology",
        "TYPE": "file",
        "EXTENSION": ("prmtop", "top", "parm7",),
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
        "CATEGORY": "structure",
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
        "CATEGORY": "trajectory",
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
        "CATEGORY": "topology",
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
        "CATEGORY": "structure",
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
        "CATEGORY": "trajectory",
        "EXTENSION": ("xtc", "trr", "nc", "ncdf", "trj", "mdcrd", "dcd"),
    },
    "PDB": {
        "REQUIRED": True,
        "TYPE": "file",
        "CATEGORY": "structure",
        "EXTENSION": ("pdb",),
    },
    "TOP": {
        "REQUIRED": False,
        "CATEGORY": "topology",
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


def get_struc_top_traj_fnames(
        system: dict, allowStructure=False,
        joinPath=None) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """Returns filenames of structure/topology/trajectory according to system's engine.

    Args:
        system (dict): Databank System
        allowStructure (bool, optional): Allow using structure instead of topology for
            further MDAnalysis' Universe initialization. Defaults to False.
        joinPath (str,None): path to be joined with filenames

    Returns ((string, string, string)): structure filename, topology filename,
            trajectory filename

    Raises:
        ValueError: function
    """
    sft = system["SOFTWARE"]
    sftSpec = software_dict[sft.upper()]
    trjF = None
    for k, v in sftSpec.items():
        if ("CATEGORY" in v and v["CATEGORY"] == "trajectory" and
                k in system and system[k] is not None):
            trjF = system[k]
            break
    topF = None
    for k, v in sftSpec.items():
        if ("CATEGORY" in v and v["CATEGORY"] == "topology" and
                k in system and system[k] is not None):
            topF = system[k]
            break
    strF = None
    for k, v in sftSpec.items():
        if ("CATEGORY" in v and v["CATEGORY"] == "structure"
                and k in system and system[k] is not None):
            strF = system[k]
            break

    def is_sequence(var):
        return isinstance(var, Sequence) and not isinstance(var, (str, bytes))

    # taking first elements
    # TODO can it be removed?
    if is_sequence(topF):
        topF = topF[0]
        if is_sequence(topF):
            topF = topF[0]

    if is_sequence(trjF):
        trjF = trjF[0]
        if is_sequence(trjF):
            trjF = trjF[0]

    if is_sequence(strF):
        strF = strF[0]
        if is_sequence(strF):
            strF = strF[0]

    # applying joinPath
    if joinPath is not None:
        if topF is not None:
            topF = os.path.join(joinPath, topF)
        if trjF is not None:
            trjF = os.path.join(joinPath, trjF)
        if strF is not None:
            strF = os.path.join(joinPath, strF)

    return strF, topF, trjF
