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

import gc
import os, sys
import warnings

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from DatabankLib import NMLDB_SIMU_PATH
from DatabankLib.analyze_nmrpca import *
from DatabankLib.databankLibrary import initialize_databank, lipids_dict

if __name__ == "__main__":
    # TODO: do we need these warnings?
    warnings.filterwarnings("ignore", category=UserWarning)

    # searches through every subfolder of a path
    # and finds every trajectory in databank
    systems = initialize_databank()

    i = 0
    listall = []

    eq_time_fname = "eq_times.json"

    for readme in systems:
        print('== Doi: ' + readme['DOI'])
        print('== System name: ' + readme['SYSTEM'])
        # getting data from databank and preprocessing them
        # Start Parser
        # TODO: 2test|    parser = Parser(NMLDB_SIMU_PATH, readme, eq_time_fname, testTraj)
        parser = Parser(NMLDB_SIMU_PATH, readme, eq_time_fname)
        # Check trajectory
        print(NMLDB_SIMU_PATH, testTraj, readme['path'])
        print(parser.indexingPath)
        print(parser.validatePath())
        if parser.validatePath() < 0:
            continue
        # Download files

        eq_time_path = os.path.join(NMLDB_SIMU_PATH, readme['path'], "eq_times.json")
        print(eq_time_path)
        if os.path.isfile(eq_time_path):
            continue

        if ('WARNINGS' in readme and 
            readme['WARNINGS'] is not None and 
            'AMBIGUOUS_ATOMNAMES' in readme['WARNINGS']):
            continue

        if ('WARNINGS' in readme and 
            readme['WARNINGS'] is not None and
            'SKIP_EQTIMES' in readme['WARNINGS']):
            continue

        parser.downloadTraj()
        # Prepare trajectory
        parser.prepareTraj()
        # Concatenate trajectory
        parser.concatenateTraj()
        equilibration_times = {}
        # Iterate over trajectories for different lipids
        for traj in parser.concatenated_trajs:
            print(f"Main: parsing lipid {traj[4]}")
            # Creat PCA for trajectory
            pca_runner = PCA(traj[0], traj[1], traj[2], traj[3], parser.trjLen)
            # Run PCA
            data = pca_runner.PCA()
            print("Main: PCA: done")
            # Project trajectory on PC1
            pca_runner.get_proj(data)
            print("Main: Projections: done")
            # Calculate autocorrelations
            pca_runner.get_autocorrelations()
            print("Main: Autocorrelations: done")
            # Estimate equilibration time
            te2 = TimeEstimator(pca_runner.autocorrelation).calculate_time()
            equilibration_times[traj[4]] = te2 / parser.trjLen
            print("Main: EQ time: done")

            print(te2 / parser.trjLen)
        gc.collect()
        parser.dumpData(equilibration_times)
