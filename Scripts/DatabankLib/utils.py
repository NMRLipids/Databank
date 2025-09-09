"""
Utility functions for Datbank standalone scripts.

Currently only ``run_analysis``, which runs a given analysis method on a
range of systems in the databank.
"""
from collections.abc import Callable
from logging import Logger

from DatabankLib import RCODE_COMPUTED, RCODE_ERROR, RCODE_SKIPPED
from DatabankLib.core import initialize_databank


def run_analysis(
    method: Callable,
    logger: Logger,
    id_range=(None, None),
) -> None:
    """
    Apply analysis ``method`` to the entire databank.

    :param method: (Callable) will be called as ``fun(system, logger)``
    :param logger: (Logger) reference to Logger initialized by the top script
    :param id_range: (A,B) filter for systems to analyze, default is
                     (None, None) which means all systems. Can be also (None, -1)
                     which means all new systems, or (0, None) which means all
                     old systems.

    :return: None
    """
    systems = initialize_databank()

    list_ids = [s["ID"] for s in systems]
    if id_range[0] is not None:
        list_ids = [s for s in list_ids if s >= id_range[0]]
    if id_range[1] is not None:
        list_ids = [s for s in list_ids if s <= id_range[1]]
    logger.info("Filtering systems by range: %s", str(id_range))
    logger.info("Number of systems in databank: %s", str(len(list_ids)))

    result_dict = {RCODE_COMPUTED: 0, RCODE_SKIPPED: 0, RCODE_ERROR: 0}

    for sid in list_ids:
        system = systems.loc(sid)
        logger.info("System title: %s", system["SYSTEM"])
        logger.info("System path: %s", system["path"])
        res = method(system, logger)
        result_dict[res] += 1

    print(f"""
    COMPUTED: {result_dict[RCODE_COMPUTED]}
    SKIPPED: {result_dict[RCODE_SKIPPED]}
    ERROR: {result_dict[RCODE_ERROR]}
    """)
