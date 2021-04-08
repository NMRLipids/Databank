# This is the NMRlipids databank containing MD simulations lipid bilayers created by the [NMRlipids project](http://nmrlipids.blogspot.com/2021/03/second-online-meeting-on-nmrlipids.html). 

This is overlay databank containing all essential information on simulations including the permanent location of simulation files, but not the trajectories themself. The databank format enables completely automatic analysis of contributed data. 

Codes used to [add](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank) and [analyze](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank) data can be found from the [scripts folder](https://github.com/NMRLipids/Databank/tree/main/Scripts). The data entries are in [Data/simulations](https://github.com/NMRLipids/Databank/tree/main/Data/Simulations) folder under subfolders named based on file hash identities. 

## Instructions to add data

1. Clone this repo to your own computer.
2. Make a new directory with the next free integer into [Scripts/BuildDatabank/info_files/](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank/info_files) folder.
3. Create info.yaml file into the folder that you created following the instructions.
4. Return to the [Databank/Scripts/BuildDatabank/](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank) folder and run
`python3 AddData.py -f {path to the info file that you created}`.
After this is finished, you should see a new folder in [Data/simulations](https://github.com/NMRLipids/Databank/tree/main/Data/Simulations) which contains the README.yaml file of your system and calculated order parameters.
5. Commit the created README.yaml and order parameter files, and make a pull request to the master branch.

## Instructions to analyze data

Simulations in the databank can be browsed going through all the files in [Data/simulations](https://github.com/NMRLipids/Databank/tree/main/Data/Simulations) folder. Jupyter-notebook template and other examples can be found from [Scripts/AnalyzeDatabank](https://github.com/NMRLipids/Databank/tree/main/Scripts/AnalyzeDatabank) folder.
