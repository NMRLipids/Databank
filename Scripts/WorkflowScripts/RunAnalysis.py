from DatabankLib import NMLDB_ROOT_PATH
from WorkflowScripts.Workflow_utils import *  
"""
Executes pipeline of analysis methods that work globally on simulations.

.. note::
   This file is only meant to be used by automated workflows.
   Users of the Databank repository can safely ignore it.
"""

def main():    
    path_dict = get_databank_paths(NMLDB_ROOT_PATH)
    run_python_script(path_dict["searchDATABANK_path"], error_message="SearchDatabank failed")
    run_python_script(path_dict["QualityEvaluation_path"], error_message="QualityEvaluation failed")
    run_python_script(path_dict["makeRanking_path"], error_message="makeRanking failed")

if __name__ == "__main__":
    main()