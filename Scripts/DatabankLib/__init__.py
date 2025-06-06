"""
For now we just need to declare this folder to be package.
TODO: organize proper package structure
"""

# TODO: think which modules should be imported at package-level
# from . import some_module

import os

NMLDB_ROOT_PATH: str = os.environ.get(
    "NMLDB_ROOT_PATH",
    os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
)
""" Path to the project root """

NMLDB_DATA_PATH: str = os.environ.get(
    "NMLDB_DATA_PATH",
    os.path.join(NMLDB_ROOT_PATH, 'Data')
)
""" Path to the project data folder """

NMLDB_SIMU_PATH: str = os.environ.get(
    "NMLDB_SIMU_PATH",
    os.path.join(NMLDB_DATA_PATH, 'Simulations')
)
""" Path to the project simulations folder"""

NMLDB_MOL_PATH: str = os.path.join(NMLDB_DATA_PATH, 'Molecules')
""" Path to the project molecules folder """

NMLDB_EXP_PATH: str = os.path.join(NMLDB_DATA_PATH, 'experiments')
""" Path to the project experiments folder """

if not os.path.isdir(NMLDB_DATA_PATH):
    raise RuntimeError(
            "Seems that you installed package in a non-debug mode. "
            "In this case you *must* set NMLDB_ROOT_PATH explicitly "
            "to git-clonned `Databank` folder.")

# Universal success codes

RCODE_SKIPPED: int = 0
""" Success code 0: calculation skipped """

RCODE_COMPUTED: int = 1
""" Success code 1: calculation successful """

RCODE_ERROR: int = 2
""" Success code 2: calculation failed """
