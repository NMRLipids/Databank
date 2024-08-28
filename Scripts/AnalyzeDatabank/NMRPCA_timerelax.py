#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main runs the loop over all databank trajectories and does the following:
    1. Calls the parser to check, merge, and concatenate trajectory
    2. Calls PCA to gep PCs, projections, and ACs
    3. Calls Time calculation

The parser assumes that the trajectory is downloaded and lipids are made
whole, or GROMACS is installed and available via gmx command.

For details on NMRPCA, see `DatabankLib/analyze_nmrpca.py`.
"""

import os, sys

import logging
logger = logging.getLogger(__name__)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import DatabankLib
from DatabankLib.core import initialize_databank
from DatabankLib.analyze import computeNMRPCA


systems = initialize_databank()
resDict = {DatabankLib.RCODE_COMPUTED: 0, 
           DatabankLib.RCODE_SKIPPED: 0,
           DatabankLib.RCODE_ERROR: 0}

for system in systems:
    logger.info("System title: " + system['SYSTEM'])
    logger.info("System path: " + system['path'])
    res = computeNMRPCA(system)
    resDict[res] += 1

print(f"""
COMPUTED: {resDict[DatabankLib.RCODE_COMPUTED]}
 SKIPPED: {resDict[DatabankLib.RCODE_SKIPPED]}
   ERROR: {resDict[DatabankLib.RCODE_ERROR]}
""")