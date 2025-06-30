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
    work_directory = "/tmp/databank_workdir"
    os.makedirs(work_directory, exist_ok=True)

    run_python_script(
        path_dict["AddData_path"],
        args=["-f", info_file_path, "-w", work_directory],
        error_message="AddData failed"
    )
    run_command(
        path_dict["calcProperties_path"],
        "Calcproperties failed",
        working_dir=path_dict["AnalyzeDatabank_path"]
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