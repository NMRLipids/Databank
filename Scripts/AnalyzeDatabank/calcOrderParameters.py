#!/usr/bin/env python3
# coding: utf-8

from DatabankLib.utils import run_analysis
from DatabankLib.analyze import computeOP
import logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    run_analysis(computeOP, logger)
