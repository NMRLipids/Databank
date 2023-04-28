# NMRlipids Databank 
NMRlipids databank is a community-driven catalogue containing atomistic MD simulations of biologically relevant lipid membranes emerging from the [NMRlipids open collaboration](http://nmrlipids.blogspot.com/2021/03/second-online-meeting-on-nmrlipids.html). 

NMRlipids databank is an overlay databank. For each simulation, there is a README.yaml which contains all the essential information for the data upcycling and reuse, including the permanent location of each simulation file. The description of the content in README.yaml files can be found from [here](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/info_files/README.md). The README.yaml files are located in [Data/simulations](https://github.com/NMRLipids/Databank/tree/main/Data/Simulations) folder under subfolders named based on file hash identities. 

## NMRlipids Databank-GUI
NMRlipids Databank-GUI, available at databank.nmrlipids.fi, provides easy access to the NMRlipids Databank content
through a graphical user interface (GUI). Simulations can be searched based on their molecular composition, force field,
temperature, membrane properties, and quality; the search results are ranked based on the simulation quality as evaluated
against experimental data when available. Membranes can be visualized, and properties between different simulations and
experiments compared.

## NMRlipids Databank-API
The databank format enables completely automatic analysis of contributed data. 

Codes used to [add](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank) and [analyze](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank) data can be found from the [scripts folder](https://github.com/NMRLipids/Databank/tree/main/Scripts). 

## Instructions to add data

1. Clone this repo to your own computer.
2. Make a new directory with the next free integer into [Scripts/BuildDatabank/info_files/](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank/info_files) folder.
3. Create info.yaml file into the folder with the help of [instructions](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/info_files/README.md).
4. Return to the [Databank/Scripts/BuildDatabank/](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank) folder and run
`python3 AddData.py -f {path to the info file that you created}`.
After this is finished, you should see a new folder in [Data/simulations](https://github.com/NMRLipids/Databank/tree/main/Data/Simulations) which contains the README.yaml file of your system and calculated order parameters.
5. Commit the created README.yaml and order parameter files, and make a pull request to the master branch.

## Instructions to analyze data

Simulations in the databank can be browsed going through all the files in [Data/simulations](https://github.com/NMRLipids/Databank/tree/main/Data/Simulations) folder. Jupyter-notebook template and other examples can be found from [Scripts/AnalyzeDatabank](https://github.com/NMRLipids/Databank/tree/main/Scripts/AnalyzeDatabank) folder.
