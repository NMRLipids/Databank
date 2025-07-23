from DatabankLib import NMLDB_ROOT_PATH
from WorkflowScripts.Workflow_utils import *  
import os 
import argparse
"""
Executes standard pipeline of processing for an info file.

.. note::
   This file is only meant to be used by automated workflows.
   Users of the Databank repository can safely ignore it.
"""

def main(info_file_path): 
    path_dict = get_databank_paths(NMLDB_ROOT_PATH)
    parent_folder = os.path.dirname(NMLDB_ROOT_PATH)

    base_tmp = os.path.join(parent_folder, "databank_workdir")
    work_directory_dry  = os.path.join(base_tmp, "dry")
    work_directory_real = os.path.join(base_tmp, "real")
    try:
        os.makedirs(work_directory_dry, exist_ok=True)
        os.makedirs(work_directory_real, exist_ok=True)
    except OSError as e:
        logger.error(f"Failed to create work directories")
        sys.exit(1)

    run_python_script(
        path_dict["AddData_path"],
        args=["-f", info_file_path, "-w", work_directory_dry, "--dry-run"],
        error_message="AddData dry run failed"
    )
    run_python_script(
        path_dict["AddData_path"],
        args=["-f", info_file_path, "-w", work_directory_real],
        error_message="AddData failed"
    )
    run_python_script(
        path_dict["compute_databank_path"],
        args = ["--nmrpca", "--maicos", "--op", "--thickness","--apl", "--range", "*-0"],
        error_message="Compute_databank failed"
    )
    delete_info_file(info_file_path)

#Gets arguments from parser
def get_args():
    parser = argparse.ArgumentParser(description="Run NMRLipids data pipeline on a YAML info file.")
    parser.add_argument("--info_file_path", required=True, help="Path to the info.yml file")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    info_file_path = args.info_file_path
    if not info_file_path:
        print("No path provided, exiting.")
        sys.exit(0)
    main(info_file_path)