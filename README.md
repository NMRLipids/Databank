# NMRlipids Databank 
This is the NMRlipids databank &mdash; a community-driven catalogue containing atomistic MD simulations of biologically relevant lipid membranes emerging from the [NMRlipids open collaboration](http://nmrlipids.blogspot.com/2021/03/second-online-meeting-on-nmrlipids.html). 

NMRlipids databank is an overlay databank. For each simulation, there is a README.yaml which contains all the essential information for the data upcycling and reuse, including the permanent location of each simulation file. The description of the content in README.yaml files can be found from [here](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/info_files/README.md). The README.yaml files are located in [Data/simulations](https://github.com/NMRLipids/Databank/tree/main/Data/Simulations) folder under subfolders named based on file hash identities. 

More information is available from the [NMRlipids databank manuscript](https://doi.org/10.26434/chemrxiv-2023-jrpwm).

## NMRlipids Databank-GUI
[NMRlipids Databank-GUI](https://databank.nmrlipids.fi/) provides easy access to the NMRlipids Databank content
through a graphical user interface (GUI). Simulations can be searched based on their molecular composition, force field,
temperature, membrane properties, and quality; the search results are ranked based on the simulation quality as evaluated
against experimental data when available. Membranes can be visualized, and properties between different simulations and
experiments compared.

## NMRlipids Databank-API
[The NMRlipids Databank-API](https://github.com/NMRLipids/Databank/tree/main/Scripts) provides programmatic access to all simulation data in the NMRlipids Databank. This enables wide range of novel data-driven applications &mdash; from construction of machine learning models that predict membrane properties, to automatic analysis of virtually any property across all simulations in the Databank. 

Minimum example for looping over available simulations is available in [here](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/template.ipynb). Codes that analyze [area per lipid](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/calcAPL.py), [C-H bond order parameters](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/calcOrderParameters.py), [X-ray scattering form factors](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/calc_FormFactors.py), and [principal component equilibration](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/NMRPCA_timerelax.py) are also available as examples. [Molecule](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank/info_files#composition-compulsory) and [atom naming](https://nmrlipids.blogspot.com/2022/04/new-yaml-format-of-mapping-files.html) conventions in each simulation are connected to the universal molecule and atom names using [mapping files](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank/mapping_files). The results can be stored using the same structure as used for README.yaml files as done, for example, for [water permeation](https://github.com/NMRLipids/DataBankManuscript/tree/main/Data/MD-PERMEATION) and lipid [flip-flop](https://github.com/NMRLipids/DataBankManuscript/tree/main/Data/Flipflops) rates in the [repository related to the NMRlipids databank manuscript](https://github.com/NMRLipids/DataBankManuscript).

## Instructions to add data into the databank

The NMRlipids Databank is open for additions of simulation data by anyone.
1. Clone this repo to your own computer.
2. Make a new directory with the next free integer into [Scripts/BuildDatabank/info_files/](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank/info_files) folder.
3. Create info.yaml file into the folder with the help of [instructions](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/info_files/README.md).
4. Return to the [Databank/Scripts/BuildDatabank/](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank) folder and run
`python3 AddData.py -f {path to the info file that you created}`.

Optional: See `python3 AddData.py --help` for available command line arguments:
```
usage: AddData.py Script [-h] [-f FILE] [-d] [-n] [-w WORK_DIR]
                         [-o OUTPUT_DIR]

Add a new dataset to the NMRLipids databank

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Input config file in yaml format.
  -d, --debug           enable debug logging output
  -n, --no-cache        always redownload repository files
  -w WORK_DIR, --work-dir WORK_DIR
                        set custom temporary working directory
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        set custom output directory
```

After this is finished, you should see a new folder in [Data/simulations](https://github.com/NMRLipids/Databank/tree/main/Data/Simulations) which contains the **README.yaml** file of your system and calculated order parameters.

5. Commit the created README.yaml file, and make a pull request to the master branch.

## System requirements

The code has been tested in Linux environment with python 3.7 or newer and recent [Gromacs](https://manual.gromacs.org/current/install-guide/index.html) version installed.

Setup using conda as distribution:

    $ conda create --name databank python==3.7.16 MDAnalysis MDAnalysisTests
    $ conda activate databank
    $ (databank) pip install tqdm pyyaml
