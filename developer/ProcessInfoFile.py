"""
Executes standard pipeline of processing for an info file.

.. note::
   This file is only meant to be used by automated workflows.
   Users of the Databank repository can safely ignore it.
"""

import argparse
import os
import sys

from WorkflowScripts.Workflow_utils import get_databank_paths, logger, run_python_script

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def run_analysis(info_file_path: str) -> None:
    """Run full analysis of info-file"""
    work_directory_real, _ = setup_folders()
    path_dict = get_databank_paths(REPO_ROOT)
    run_python_script(
        path_dict["adddata_path"],
        args=["-f", info_file_path, "-w", work_directory_real],
        error_message="AddData failed",
    )
    run_python_script(
        path_dict["compute_databank_path"],
        args=["--nmrpca", "--ff", "--op", "--thickness", "--apl", "--range", "*-0"],
        error_message="Compute_databank failed",
    )
    delete_info_file(info_file_path)


def run_dry_run(info_file_path: str) -> None:
    """Run AddData dry-run for pre-analysis of info file to rule out errors"""
    _, work_directory_dry = setup_folders()
    path_dict = get_databank_paths(REPO_ROOT)
    run_python_script(
        path_dict["adddata_path"],
        args=["-f", info_file_path, "-w", work_directory_dry, "--dry-run"],
        error_message="AddData dry run failed",
    )


def setup_folders() -> tuple[str, str]:
    """Set up folders for processing info file"""
    parent_folder = os.path.dirname(REPO_ROOT)
    base_tmp = os.path.join(parent_folder, "databank_workdir")
    work_directory_dry = os.path.join(base_tmp, "dry")
    work_directory_real = os.path.join(base_tmp, "real")
    try:
        os.makedirs(work_directory_dry, exist_ok=True)
        os.makedirs(work_directory_real, exist_ok=True)
    except OSError:
        logger.exception("Failed to create work directories")
        sys.exit(1)
    return work_directory_real, work_directory_dry


def get_args() -> argparse.Namespace:
    """Retrieve arguments from parser"""
    parser = argparse.ArgumentParser(description="Run NMRLipids data pipeline on a YAML info file.")
    parser.add_argument("--info_file_path", required=True, help="Path to the info.yml file")
    parser.add_argument("--dry-run", action="store_true", help="Only run preanalysis (no compute/commit)")
    return parser.parse_args()


def delete_info_file(info_file_path: str) -> None:
    """
    Delete an info file at the specified path.

    :param info_file_path: Path to the info file to be removed.
    :raises OSError: Prints a warning if the file cannot be deleted.
    """
    try:
        os.remove(info_file_path)
    except FileNotFoundError:
        logger.exception(f"Info file not found, nothing to delete: {info_file_path}")
    except OSError:
        logger.exception(f"Could not delete {info_file_path}")
    else:
        logger.info(f"Deleted info file: {info_file_path}")


if __name__ == "__main__":
    args = get_args()
    info_file_path = args.info_file_path
    if not info_file_path:
        logger.info("No path provided, exiting.")
        sys.exit(0)
    if args.dry_run:
        run_dry_run(info_file_path)
    else:
        run_analysis(info_file_path)
