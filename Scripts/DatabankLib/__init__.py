"""
For now we just need to declare this folder to be package.
TODO: organize proper package structure
"""

#TODO: think which modules should be imported at package-level
# from . import some_module

import os

NMLDB_ROOT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

NMLDB_DATA_PATH = os.path.join(NMLDB_ROOT_PATH, 'Data')
NMLDB_SIMU_PATH = os.path.join(NMLDB_ROOT_PATH, 'Data', 'Simulations')
NMLDB_EXP_PATH  = os.path.join(NMLDB_ROOT_PATH, 'Data', 'experiments')