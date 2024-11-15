from DatabankLib import RCODE_COMPUTED, RCODE_ERROR, RCODE_SKIPPED
from DatabankLib.core import initialize_databank


from logging import Logger
from typing import Callable


def run_analysis(method: Callable, logger: Logger):
    """Apply analysis ``method`` to the entire databank.

    Args:
        method (Callable): will be called as ``fun(system, logger)``
        logger (Logger): reference to Logger initialized by the top script
    """
    systems = initialize_databank()
    resDict = {
        RCODE_COMPUTED: 0,
        RCODE_SKIPPED: 0,
        RCODE_ERROR: 0}

    for system in systems:
        logger.info("System title: " + system['SYSTEM'])
        logger.info("System path: " + system['path'])
        res = method(system, logger)
        resDict[res] += 1

    print(f"""
    COMPUTED: {resDict[RCODE_COMPUTED]}
    SKIPPED: {resDict[RCODE_SKIPPED]}
    ERROR: {resDict[RCODE_ERROR]}
    """)
