"""
:module: settings/molecules.py
:description: Module file with definition of different global-level dictionaries.

There is a dictionary of lipids, ions, etc. If you add a lipid which is not yet
in the databank, you have to add it here!
"""

# Dictionary of possible lipids
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
    "TMCL": {
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
    "BOG": {
        "REQUIRED": False,
        "TYPE": "string",
    }
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
    "TMA": {
        "REQUIRED": False,
        "TYPE": "string",
    }
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
