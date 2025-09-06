"""
Execute pipeline of analysis methods that work globally on simulations.

.. note::
   This file is only meant to be used by automated workflows.
   Users of the Databank repository can safely ignore it.
"""

from DatabankLib import NMLDB_ROOT_PATH
from WorkflowScripts.Workflow_utils import get_databank_paths, run_python_script


def main() -> None:
    """Run analysis scripts on simulations."""
    path_dict = get_databank_paths(NMLDB_ROOT_PATH)
    run_python_script(path_dict["searchdatabank_path"], error_message="SearchDatabank failed")
    run_python_script(path_dict["qualityevaluation_path"], error_message="QualityEvaluation failed")
    run_python_script(path_dict["makeranking_path"], error_message="MakeRanking failed")


if __name__ == "__main__":
    main()
