"""
Contains methods used for python scripts related to workflows.

.. note::
   This module is only used by automated workflows. Users of the Databank
   repository can safely ignore it.
"""

import logging
import os
import subprocess
import sys

logging.basicConfig(
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%I:%M:%S %p",
    level=logging.INFO,
)
logger = logging.getLogger()


def run_python_script(
    script_path: str,
    args: list | None = None,
    error_message: str = "Python script failed",
    working_dir: str | None = None,
) -> None:
    """
    Execute a Python script with the current interpreter and optional arguments.

    :param script_path: Absolute path to the Python script to run.
    :param args: List of arguments to pass to the script (defaults to []).
    :param error_message: Message to display if execution fails.
    :param working_dir: Optional working directory in which to run the script.
    :raises SystemExit: Exits with code 1 if the script execution fails.
    """
    if args is None:
        args = []

    if not os.path.isfile(script_path):
        logger.error(f"Script not found: {script_path}")
        sys.exit(1)

    try:
        logger.info(f"Running python script with path {script_path}")
        subprocess.run(
            [sys.executable, script_path, *args],
            check=True,
            cwd=working_dir,
        )
    except OSError:
        logger.exception(f"{error_message}, caught OSError")
        sys.exit(1)
    except ValueError:
        logger.exception(f"{error_message}, caught ValueError")
        sys.exit(1)
    except subprocess.SubprocessError:
        logger.exception(f"{error_message}, caught SubprocessError:")
        sys.exit(1)


def get_databank_paths(nmlb_root_path: str) -> dict:
    """
    Retrieve relevant paths from databank.

    :param nmlb_root_path: Root path of the database repository.
    :return: Dictionary mapping descriptive keys to full paths of databank components.
    :rtype: dict
    """
    builddatabank_path = os.path.join(nmlb_root_path, "Scripts", "BuildDatabank")
    adddata_path = os.path.join(builddatabank_path, "AddData.py")
    analyzedatabank_path = os.path.join(nmlb_root_path, "Scripts", "AnalyzeDatabank")
    compute_databank_path = os.path.join(analyzedatabank_path, "compute_databank.py")
    searchdatabank_path = os.path.join(builddatabank_path, "searchDATABANK.py")
    qualityevaluation_path = os.path.join(builddatabank_path, "QualityEvaluation.py")
    makeranking_path = os.path.join(builddatabank_path, "makeRanking.py")

    return {
        "builddatabank_path": builddatabank_path,
        "adddata_path": adddata_path,
        "analyzedatabank_path": analyzedatabank_path,
        "compute_databank_path": compute_databank_path,
        "searchdatabank_path": searchdatabank_path,
        "qualityevaluation_path": qualityevaluation_path,
        "makeranking_path": makeranking_path,
    }


