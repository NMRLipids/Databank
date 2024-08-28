#!/usr/bin/env python3
# coding: utf-8

import os, sys

import logging
logger = logging.getLogger(__name__)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import DatabankLib
from DatabankLib.core import initialize_databank
from DatabankLib.analyze import computeOP


systems = initialize_databank()
resDict = {DatabankLib.RCODE_COMPUTED: 0, 
           DatabankLib.RCODE_SKIPPED: 0,
           DatabankLib.RCODE_ERROR: 0}

for system in systems:
    logger.info("System title: " + system['SYSTEM'])
    logger.info("System path: " + system['path'])
    res = computeOP(system)
    resDict[res] += 1

print(f"""
COMPUTED: {resDict[DatabankLib.RCODE_COMPUTED]}
 SKIPPED: {resDict[DatabankLib.RCODE_SKIPPED]}
   ERROR: {resDict[DatabankLib.RCODE_ERROR]}
""")
