# NMRlipids Databank 

This is the NMRlipids databank &mdash; a community-driven catalogue containing atomistic MD simulations of biologically relevant lipid membranes emerging from the [NMRlipids open collaboration](http://nmrlipids.blogspot.com/2021/03/second-online-meeting-on-nmrlipids.html). 


# Documentation

The NMRlipids databank documentation is available in [here](https://nmrlipids.github.io/index.html). 
More information and example applications are available from the [NMRlipids databank manuscript](https://doi.org/10.1038/s41467-024-45189-z).

## API

The `DatabankLib` python module provides programmatic access to all simulation data in the 
NMRlipids Databank. It allows to request data to construct various datasets on the base of 
experimental and simulation data from the Databank that allow one to learn various models about 
lipid bilayer properties. It also allows to design and effectively run automated analysis across
all the simulations in the Databank.

The documentation for `DatabankLib` module is available in [here](https://nmrlipids.github.io/auto_gen/Scripts.DatabankLib.html).

## How to use 

A [jupyter template notebook](https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb) can be used to get started with the analyses utilizing the NMRlipids databank.

Connection of [Universal molecule and atom naming conventions](https://nmrlipids.github.io/moleculesAndMapping.html) with simulation specific names delivered by mapping files can be used to perform automatic analyses over large sets of simulations. The results for large analyses can be stored using the same structure as used for `README.yaml` files as done, for example, for [water permeation](https://github.com/NMRLipids/DataBankManuscript/tree/main/Data/MD-PERMEATION) and lipid [flip-flop](https://github.com/NMRLipids/DataBankManuscript/tree/main/Data/Flipflops) rates in the [repository related to the NMRlipids databank manuscript](https://github.com/NMRLipids/DataBankManuscript).


# GUI

[NMRlipids Databank-GUI](https://databank.nmrlipids.fi/) provides easy access to the NMRlipids Databank content
through a graphical user interface (GUI). 
Simulations can be searched based on their molecular composition, force field,
temperature, membrane properties, and quality; the search results are ranked based on the simulation quality as evaluated
against experimental data when available. GUI provides basic graphical reports
for the computed properties as well as graphical comparison between simulation
and experimental data.

The GUI is being developed in the repository [BilayerGUI_laravel](https://github.com/NMRlipids/BilayerGUI_laravel).


# Installation

The code has been tested in Linux and MacOS environment with python 3.10 or newer. Recent [Gromacs](https://manual.gromacs.org/current/install-guide/index.html) version should be available in the system. All dependecies are listed in [requirements.txt](Scripts/DatabankLib/requirements.txt).

Note that the data is stored as a submodule repository and should be loaded after clonning. Default data storage is [BilayerData](https://github.com/NMRLipids/BilayerData), and it is loaded automatically by using
```
$ git submodule update --init --remote
```

We recomend installing python libraries into an environment, for example, using conda:

```
 $ conda create --name databank python==3.10 MDAnalysis periodictable -c conda-forge
 $ conda activate databank
 $ (databank) conda install tqdm yaml -c conda-forge
```

You should also activate *DatabankLib* package:

```
 $ (databank) cd Databank
 $ (databank) pip install -e .
```

You can install the package in non-development mode, without `-e` (it's obligatory in Colab runtime environment); however, in this case, the package will be installed to the folder with other pip-packages and it will not know about the path to the `Data` folder. Then you should provide the path to the *repository root* by setting the environment variable `NMLDB_ROOT_PATH`.


# Contribution

The project is open for contributions! 

For code development, please use extended requirements described in `requirements-dev.txt`:
```
 $ (databank) pip install -e . -r Scripts/DatabankLib/requirements-dev.txt
```
It will install `pytest` for unit tests and `ruff` for syntax check.

Please consult [CONTRIBUTION.md](./CONTRIBUTION.md) for further information.