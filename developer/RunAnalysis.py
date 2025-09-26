"""
Execute pipeline of analysis methods that work globally on simulations.

.. note::
   This file is only meant to be used by automated workflows.
   Users of the Databank repository can safely ignore it.
"""

from WorkflowScripts.Workflow_utils import get_databank_paths, run_python_script

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def main() -> None:
    """Run analysis scripts on simulations."""
    path_dict = get_databank_paths(REPO_ROOT)
    run_python_script(path_dict["searchdatabank_path"], error_message="SearchDatabank failed")
    run_python_script(path_dict["qualityevaluation_path"], error_message="QualityEvaluation failed")
    run_python_script(path_dict["makeranking_path"], error_message="MakeRanking failed")


if __name__ == "__main__":
    main()
