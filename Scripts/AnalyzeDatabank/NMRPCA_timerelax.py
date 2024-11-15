#!/usr/bin/env python3
# coding: utf-8
"""
:program: NMRPCA_timerelax.py
:description: Compute PCA-derived lipid relaxation times.

Main runs the loop over all databank trajectories and does the following:
    1. Calls the parser to check, merge, and concatenate trajectory
    2. Calls PCA to gep PCs, projections, and ACs
    3. Calls Time calculation

The parser assumes that the trajectory is downloaded and lipids are made
whole, or GROMACS is installed and available via gmx command.

For details on NMRPCA, see `DatabankLib/analyze_nmrpca.py`.
"""

from DatabankLib.utils import run_analysis
from DatabankLib.analyze import computeNMRPCA
import logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    run_analysis(computeNMRPCA, logger)
