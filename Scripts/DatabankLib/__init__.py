"""
For now we just need to declare this folder to be package.
TODO: organize proper package structure
"""

# TODO: think which modules should be imported at package-level
# from . import some_module

import os

NMLDB_ROOT_PATH = os.environ.get(
    "NMLDB_ROOT_PATH",
    os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    )

NMLDB_DATA_PATH = os.path.join(NMLDB_ROOT_PATH, 'Data')
NMLDB_SIMU_PATH = os.path.join(NMLDB_ROOT_PATH, 'Data', 'Simulations')
NMLDB_EXP_PATH = os.path.join(NMLDB_ROOT_PATH, 'Data', 'experiments')

if not os.path.isdir(NMLDB_DATA_PATH):
    raise RuntimeError(
            "Seems that you installed package in a non-debug mode. "
            "In this case you *must* set NMLDB_ROOT_PATH explicitly "
            "to git-clonned `Databank` folder.")

# Universal success codes

RCODE_SKIPPED: int = 0
""" Success code 0: calculation skipped """

RCODE_COMPUTED: int = 1
""" Success code 1: calculation sucessful """

RCODE_ERROR: int = 2
""" Success code 2: calculation failed """
