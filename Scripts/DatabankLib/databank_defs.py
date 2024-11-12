#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module file with definition of different global-level dictionaries.

There is a dictionary of lipids, ions, etc.

If you add a lipid which is not yet in the databank, you have to add it here!
"""

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
    "SAPI24": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SAPI25": {
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
    "CER180": {
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
    "SLiPC": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SM16": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "SM18": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "TOCL": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "TLCL_0H": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "GM1": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "DPPGK": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "GB3": {
        "REQUIRED": False,
        "TYPE": "string",
    },
}


# Dictionary of other than lipid molecules.
#
# If you add other than a lipid molecule which is not yet in the databank, you have to
# add it here

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
    "FFSAPI24": {
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
    "FFCER180": {
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
    "FFSM16": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFSM18": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFTOCL": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFGM1": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFDPPGK": {
        "REQUIRED": False,
        "TYPE": "string",
    },
    "FFGB3": {
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
    "EDR": {
        "REQUIRED": False,
        "TYPE": "file",
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
