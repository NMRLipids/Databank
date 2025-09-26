"""
DatabankLib package for NMRlipids Databank project.

Here we define the main package variables and constants used throughout the NMRlipids
Databank project.
This includes _Package Information_, _Paths_, and _Return Codes_.
"""

import os

from ._version import __version__

# Package Information

__author__ = "NMRlipids Open Collaboration"

__author_email__ = "samuli.ollila@helsinki.fi"

__description__ = (
    "NMRlipids Databank main package. "
    "This package contains the core functionality for managing and accessing "
    "the NMRlipids Databank."
)

__license__ = "GPL-v3"

__url__ = "https://github.com/NMRlipids/Databank"

# Global Paths

NMLDB_DATA_PATH: str = os.environ.get(
    "NMLDB_DATA_PATH",
    os.path.join(os.getcwd(), "BilayerData"),
)
""" Path to the Databank Data folder """

NMLDB_SIMU_PATH: str = os.environ.get(
    "NMLDB_SIMU_PATH",
    os.path.join(NMLDB_DATA_PATH, "Simulations"),
)
""" Path to the project simulations folder"""

NMLDB_MOL_PATH: str = os.path.join(NMLDB_DATA_PATH, "Molecules")
""" Path to the project molecules folder """

NMLDB_EXP_PATH: str = os.path.join(NMLDB_DATA_PATH, "experiments")
""" Path to the project experiments folder """

# Universal return codes

RCODE_SKIPPED: int = 0
""" Success code 0: calculation skipped """

RCODE_COMPUTED: int = 1
""" Success code 1: calculation successful """

RCODE_ERROR: int = 2
""" Success code 2: calculation failed """

print(
    f"DatabankLib {__version__} by {__author__} - {__license__}",
)

if os.path.isdir(NMLDB_DATA_PATH):
    from DatabankLib.settings import molecules

    _ = len(molecules.lipids_set)
    print(
        f"NMRlipids Databank is initialized from the folder: {NMLDB_DATA_PATH}\n"
        "---------------------------------------------------------------",
    )
else:
    msg = f"""
Error: no data folder {NMLDB_DATA_PATH}.
If Data folder was not created, please create it by using
 $ nml_initialize_data.py toy
          OR
 $ nml_initialize_data.py stable
and then specify by NMLDB_DATA_PATH environment variable."""
    raise RuntimeError(msg)
