# NMRlipids Databank 
This is the NMRlipids databank &mdash; a community-driven catalogue containing atomistic MD simulations of biologically relevant lipid membranes emerging from the [NMRlipids open collaboration](http://nmrlipids.blogspot.com/2021/03/second-online-meeting-on-nmrlipids.html). 

NMRlipids databank is an overlay databank. For each simulation, there is a README.yaml which contains all the essential information for the data upcycling and reuse, including the permanent location of each simulation file. The description of the content in README.yaml files can be found from [here](https://nmrlipids.github.io/READMEcontent.html). The README.yaml files are located in [Data/simulations](https://github.com/NMRLipids/Databank/tree/main/Data/Simulations) folder under subfolders named based on file hash identities. 

The NMRlipids databank documentation is available in [here](https://nmrlipids.github.io/index.html). 

More information and example applications are available from the [NMRlipids databank manuscript](https://doi.org/10.26434/chemrxiv-2023-jrpwm).

## NMRlipids Databank-GUI
[NMRlipids Databank-GUI](https://databank.nmrlipids.fi/) provides easy access to the NMRlipids Databank content
through a graphical user interface (GUI). Simulations can be searched based on their molecular composition, force field,
temperature, membrane properties, and quality; the search results are ranked based on the simulation quality as evaluated
against experimental data when available. Membranes can be visualized, and properties between different simulations and
experiments compared.

## NMRlipids Databank-API
The NMRlipids Databank-API provides programmatic access to all simulation data in the NMRlipids Databank. This enables wide range of novel data-driven applications &mdash; from construction of machine learning models that predict membrane properties, to automatic analysis of virtually any property across all simulations in the Databank. 

NMRlipids Databank-API documentation is available in [here](https://nmrlipids.github.io/databankLibrary.html).

[A template](https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb) can be used to get started with the analyses utilizing the NMRlipids databank. Codes that analyze [area per lipid](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/calcAPL.py), [C-H bond order parameters](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/calcOrderParameters.py), [X-ray scattering form factors](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/calc_FormFactors.py), and [principal component equilibration](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/NMRPCA_timerelax.py) are also available as examples. Connection of [Universal molecule and atom naming conventions](https://nmrlipids.github.io/moleculesAndMapping.html) with simulation specific names delivered by mapping files can be used to perform automatic analyses over large sets of simulations. The results for large analyses can be stored using the same structure as used for README.yaml files as done, for example, for [water permeation](https://github.com/NMRLipids/DataBankManuscript/tree/main/Data/MD-PERMEATION) and lipid [flip-flop](https://github.com/NMRLipids/DataBankManuscript/tree/main/Data/Flipflops) rates in the [repository related to the NMRlipids databank manuscript](https://github.com/NMRLipids/DataBankManuscript).

## Instructions to add data into the databank

The NMRlipids Databank is open for additions of simulation data by anyone. Instructions to add data are available in [here](https://nmrlipids.github.io/addingData.html).

## System requirements

The code has been tested in Linux environment with python 3.7 or newer and recent [Gromacs](https://manual.gromacs.org/current/install-guide/index.html) version installed.

Setup using conda as distribution:

    $ conda create --name databank python==3.7.16 MDAnalysis MDAnalysisTests
    $ conda activate databank
    $ (databank) pip install tqdm pyyaml
