#!/usr/bin/env python3
# coding: utf-8

"""
:program: compute_databank.py
:description: Script to compute various properties for systems of the databank.

Usage:
    ./compute_databank.py [--apl] [--nmrpca] [--maicos] [--thickness] \
        [--OP] [--range 0-1000] [--debug] [-h]

This script allows you to compute properties: APL, thickness, PCA relaxation times,
MAICOS (density profiles, water orientation, dielectric profiles, X-ray form factor).

:param range: Range of system IDs to analyze, e.g. 0-1000. The interval can be
              opened at either end, e.g. *-50 or 25-*. Use *-0 to analyze all 
              just-added systems.
:param apl: Compute APL (Area Per Lipid) for all systems.
:param nmrpca: Compute NMR PCA for all systems.
:param maicos: Compute MAICOS for all systems.
:param thickness: Compute Thickness for all systems.
:param OP: Compute Order Parameter for all systems.
:param debug: Enable debug logging.
:param help: Show this help message and exit.
"""

import logging
import argparse

from DatabankLib.utils import run_analysis

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="compute_databank.py Script",
        description="Compute computable properties"
    )
    parser.add_argument(
        "--apl",
        help="Compute APL (Area Per Lipid) for all systems",
        action="store_true",
    )
    parser.add_argument(
        "--nmrpca",
        help="Compute NMR PCA for all systems",
        action="store_true",
    )
    parser.add_argument(    
        "--maicos",
        help="Compute MAICOS for all systems",
        action="store_true",
    )
    parser.add_argument(
        "--thickness",
        help="Compute Thickness for all systems",
        action="store_true",
    )
    parser.add_argument(
        "--OP",
        help="Compute Order Parameter for all systems",
        action="store_true",
    )
    parser.add_argument(
        "-r",
        "--range",
        help="Range of system IDs to analyze, e.g. 0-1000",
        type=str,
        default=None,
        metavar="ID_RANGE"
    )
    parser.add_argument(
        "-d",
        "--debug",
        help="Enable debug logging",
        action="store_true"
    )
    args = parser.parse_args()

    # configure logging
    logging_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s",
        datefmt="%I:%M:%S %p",
        level=logging_level,
    )
    logger = logging.getLogger()

    id_range = (None, None)
    if args.range:
        # range can be *-50 or 25-* or 25-38
        if args.range == "*":
            id_range = (None, None)
        elif args.range.startswith("*-"):
            try:
                id_range = (None, int(args.range[2:]))
            except ValueError:
                logger.error("Invalid ID range format. Use 'start-end' format.")
                exit(1)
        elif args.range.endswith("-*"):
            try:
                id_range = (int(args.range[:-2]), None)
            except ValueError:
                logger.error("Invalid ID range format. Use 'start-end' format.")
                exit(1)
        else:
            # assume format is start-end
            try:
                id_range = tuple(map(int, args.range.split("-")))
            except ValueError:
                logger.error("Invalid ID range format. Use 'start-end' format.")
                exit(1)

    if args.apl:
        logger.info("Computing APL (Area Per Lipid) for all systems")
        from DatabankLib.analyze import computeAPL
        run_analysis(computeAPL, logger, id_range=id_range)
    if args.nmrpca:
        logger.info("Computing NMR PCA for all systems")
        from DatabankLib.analyze import computeNMRPCA
        run_analysis(computeNMRPCA, logger, id_range=id_range)
    if args.maicos:
        logger.info("Computing MAICOS for all systems")
        from DatabankLib.analyze import computeMAICOS
        run_analysis(computeMAICOS, logger, id_range=id_range)
    if args.thickness:
        logger.info("Computing Thickness for all systems")
        from DatabankLib.analyze import computeTH
        run_analysis(computeTH, logger, id_range=id_range)
    if args.OP:
        logger.info("Computing Order Parameter for all systems")
        from DatabankLib.analyze import computeOP
        run_analysis(computeOP, logger, id_range=id_range)