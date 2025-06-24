# NMRlipids Databank 

This is the NMRlipids databank &mdash; a community-driven catalogue containing atomistic MD simulations of biologically relevant lipid membranes emerging from the [NMRlipids open collaboration](http://nmrlipids.blogspot.com/2021/03/second-online-meeting-on-nmrlipids.html). 


# Documentation

The NMRlipids databank documentation is available in [here](https://nmrlipids.github.io/index.html). 
More information and example applications are available from the [NMRlipids databank manuscript](https://doi.org/10.1038/s41467-024-45189-z).

## API

The `DatabankLib` python module provides programmatic access to all simulation data in the NMRlipids Databank. This enables wide range of novel data-driven applications &mdash; from construction of machine learning models that predict membrane properties, to automatic analysis of virtually any property across all simulations in the Databank. 

NMRlipids Databank-API documentation is available in [here](https://nmrlipids.github.io/databankLibrary.html).

## How to use 

[A template](https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb) can be used to get started with the analyses utilizing the NMRlipids databank. Codes that analyze [area per lipid](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/calcAPL.py), [C-H bond order parameters](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/calcOrderParameters.py), [X-ray scattering form factors](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/calc_FormFactors.py), and [principal component equilibration](https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/NMRPCA_timerelax.py) are also available as examples. 

Connection of [Universal molecule and atom naming conventions](https://nmrlipids.github.io/moleculesAndMapping.html) with simulation specific names delivered by mapping files can be used to perform automatic analyses over large sets of simulations. The results for large analyses can be stored using the same structure as used for `README.yaml` files as done, for example, for [water permeation](https://github.com/NMRLipids/DataBankManuscript/tree/main/Data/MD-PERMEATION) and lipid [flip-flop](https://github.com/NMRLipids/DataBankManuscript/tree/main/Data/Flipflops) rates in the [repository related to the NMRlipids databank manuscript](https://github.com/NMRLipids/DataBankManuscript).


# GUI

[NMRlipids Databank-GUI](https://databank.nmrlipids.fi/) provides easy access to the NMRlipids Databank content
through a graphical user interface (GUI). 
Simulations can be searched based on their molecular composition, force field,
temperature, membrane properties, and quality; the search results are ranked based on the simulation quality as evaluated
against experimental data when available. Membranes can be visualized, and properties between different simulations and
experiments compared.


# Installation

The code has been tested in Linux and MacOS environment with python 3.9 or newer. Recent [Gromacs](https://manual.gromacs.org/current/install-guide/index.html) version should be available in the system. All dependecies are listed in [requirements.txt](Scripts/DatabankLib/requirements.txt).

Note that the data is stored as a submodule repository and should be loaded after clonning. Default data storage is [BilayerData](https://github.com/NMRLipids/BilayerData), and it is loaded automatically by using
```
$ git submodule update --init --remote
```

We recomend installing python libraries into an environment, for example, using conda:

```
 $ conda create --name databank python==3.10 'numpy<2.0' MDAnalysis periodictable -c conda-forge
 $ conda activate databank
 $ (databank) conda install tqdm yaml -c conda-forge
```

You should also activate *DatabankLib* package:

```
 $ (databank) cd Databank
 $ (databank) pip install -e .
```

You can install the package in non-development mode, without `-e` (it's obligatory in Colab runtime environment); however, in this case, the package will be installed to the folder with other pip-packages and it will not know about the path to the `Data` folder. Then you should provide the path to the *repository root* by setting the environment variable `NMLDB_ROOT_PATH`.

# Installation using Docker

We provide a Docker-based development environment that allows for easy testing of new code and features. This method is recommended for development and testing purposes.

## 1. Install Docker

Before using the Docker-based development environment, you'll need to have Docker installed on your system:

- **Linux**: Follow the official Docker installation guide for Linux: [https://docs.docker.com/engine/install/](https://docs.docker.com/engine/install/)
- **macOS**: Download and install Docker Desktop for Mac from [https://www.docker.com/products/docker-desktop/](https://docs.docker.com/desktop/setup/install/mac-install/)
- **Windows**: Download and install Docker Desktop for Windows from [https://www.docker.com/products/docker-desktop/](https://docs.docker.com/desktop/setup/install/windows-install/) (requires Windows Subsystem for Linux)

## 2. Using the Docker Development Environment

### Setting Up the Development Environment
Go to the directory where you have the Databank repository. Then, 

1. Download the latest nmrlipids core image:
   ```
   docker pull nmrlipids/core:latest
   ```

2. Alternatively the Docker image can be built locally:
   ```bash
   docker build -t NAME_OF_THE_DOCKER_IMAGE .
   ```


2. Initialize the submodule data:
   ```bash
   git submodule update --init
   ```

### Testing Code

1. Start the container with your code mounted:
   ```bash
   docker run -it -v $(pwd):/workspace NAME_OF_THE_DOCKER_IMAGE
   ```
   or
   if you pulled the latest image,
   ```bash
   docker run -it -v $(pwd):/workspace nmrlipids/core:latest
   ```

2. Inside the container:
   ```bash
   # Install base requirements and the DatabankLib
   pip install -e . -r Scripts/DatabankLib/requirements.txt
   
   # Run tests
   ./runtests.sh
   ```
   Do not forget to add this new package to the requirements.txt! 

3. When done testing:
   ```bash
   # Exit container
   exit
   ```

### Managing Different Test Environments

By default the core image features the dev requirements from this repository but you can easily create different testing environments with different dependencies:

1. Start container again:
   ```bash
   docker run -it -v $(pwd):/workspace NAME_OF_THE_DOCKER_IMAGE
   ```

2. Inside container, create a new virtual environment using Conda as previously described. 

### Tips for Docker Development

- Your local code changes are automatically reflected in the container since we use a volume mount
- You can run multiple test environments simultaneously by creating different virtual environments
- To see all running containers: `docker ps`
- To stop a container: `docker stop <container_id>`
- To remove all resources related to docker including images, container, build-resources: ` docker system prune -a`

# Contribution

The project is open for contributions! 
For code development, please use extended requirements described in `requirements-dev.txt`:
```
 $ (databank) pip install -e . -r Scripts/DatabankLib/requirements-dev.txt
```
It will install `pytest` for unit tests and `flake8` for syntax check.
