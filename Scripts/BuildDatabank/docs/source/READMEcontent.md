(readmecontent)=

# User input and content of README.yaml files

Each simulation in the NMRlipids databank is assigned with a README.yaml file which contains all the essential information of the simulation. These files are created from the manually contributed info.yaml files as described in [Adding simulations into the NMRlipids databank](addData). The README.yaml files are stored in the [NMRlipids GitHub repository](https://github.com/NMRLipids/Databank/tree/main/Data/Simulations). 

README.yaml files contain information that is manually entered into info.yaml files and automatically exctracted information by the ``AddData.py`` program. Table below lists the manually entered compulsory and optional parameters, as well as automatically extracted information from simulation files.  

key | description | type  
----|--------|------
DOI | DOI from where the raw data is found | user given (compulsory) 
TRJ | Name of the trajectory file found from DOI | 
TPR | Name of the topology file found from DOI | 
SOFTWARE | Software used to run the simulation | 
PREEQTIME | Pre-equilibrate time in nanoseconds. | 
TIMELEFTOUT | Equilibration period in the uploaded trajectory. |
COMPOSITION | Molecules names and mapping files. | 
DIR\_WRK | Temporary local working directory | 
UNITEDATOM\_DICT | Hydrogen information for united atom simulations | 
TYPEOFSYSTEM | Lipid bilayer or something else | 
||
PUBLICATION | Reference to a publication(s) related to the data. | User given (optional)
AUTHORS\_CONTACT | Name and email of the main author(s) of the data. | 
SYSTEM | System description in the free text format | 
SOFTWARE\_VERSION | Version of the used software | 
FF | Name of the used force field | 
FF\_SOURCE | Source of the force field parameters | 
FF\_DATE |  Date when force field parameters were accessed | 
FF{molename} | Molecule specific force field information. | 
CPT | Name of the Gromacs checkpoint file. | 
LOG | Name of the Gromacs log file. | 
TOP | Name of the Gromacs top file. | 
GRO | Name of the Gromacs gro file. | 
||
TRAJECTORY\_SIZE | Size of the trajectory file in bytes | automatically extracted data. 
TRJLENGTH | Lenght of the trajectory (ps). | 
TEMPERATURE | Temperature of the simulation. | 
NUMBER\_OF\_ATOMS | Number of atoms in the simulation. | 
DATEOFRUNNIG | Date when added into the databank | 
EXPERIMENT | Potentially connected experimental data | 
COMPOSITION | Numbers of molecules in the system | 

#### DOI (compulsory)
Give the DOI identity for the location of simulation files. 
Current databank works only for the data in [Zenodo](www.zenodo.org), but other potential sources are may be implemented in the future. 
Note that the DOI must point to a specific version of dataset in Zenodo. DOIs pointing to all versions of certain dataset do not work.

#### TRJ (compulsory)
Give the name of the trajectory file that is found from the DOI given above.

#### TPR (compulsory)
Give the name of the file with topology information (tpr file in the case of Gromacs) that is found from the DOI given above.

#### SOFTWARE (compulsory)
Give the name of software used to run the simulation. The options are GROMACS, AMBER, NAMD, CHARMM and OPENMM. So far, only simulations run with GROMACS are accepted by the script.


#### PREEQTIME (compulsory)
Give the time simulated before the uploaded trajectory in nanoseconds. For example, if you upload 100-200 ns part of total 200 ns simulation, this should value should be 100.

#### TIMELEFTOUT (compulsory)
Give the time that should be considered as an equilibration period in the uploaded trajectory. Frames before the give time will be discarded in the analysis. 
For example, if you upload 0-200 ns part of total 200 ns simulation where the first 100 ns should be considered as an equilibration, this value should be 100.

#### COMPOSITION (compulsory)
Information about the composition (i.e., the number of molecules) of each simulation is stored in COMPOSITION as python dictionary format.
As an input, the COMPOSITION requires the universal name for each molecule present in the simulation (for definitions, see [Universal molecule and atom names](molecule_names)) as the first keys.
For each molecule, another dictionary containing the simulation specific residue name (NAME) and the name of the mapping file (MAPPING) needs to be then defined.
Mapping file defines the universal atom names for each molecule, for details see [Universal molecule and atom names](molecule_names)).
For example, the COMPOSITION dictionary input for [a system](https://doi.org/10.5281/zenodo.259392) containing POPC, Cholesterol (CHOL), water (SOL), sodium (SOD), and chloride (CLA) is given as:

    COMPOSITION:
     CHOL:                                      # universal molecule name
      NAME: CHL1                                # simulation specific molecule name
      MAPPING: mappingCHOLESTEROLcharmm.txt     # name of the mapping file
     POPC:
      NAME: POPC
      MAPPING: mappingPOPCcharmm.txt
     SOL:
      NAME: TIP3
      MAPPING: mappingTIP3PCHARMMgui.txt
     SOD:
      NAME: SOD
      MAPPING: mappingSOD.txt
     CLA:
      NAME: CLA
      MAPPING: mappingCLA.txt

When running ``AddData.py``, the numbers of molecules are added in additional ``COUNT`` keys into dictionaries of each molecule, see COMPOSITION (ouput) below. 

#### DIR\_WRK (compulsory)
Give the path of the working directory in your local computer. The trajectory and topology files will be downloaded to this trajectory, and temporary files created during processing will be stored here. 
#### UNITEDATOM\_DICT (compulsory for united atom trajectories)
Order parameters from united atom simulations are calculated using the [buildH code](https://github.com/patrickfuchs/buildH). 
For united atom simulations, you need to tell how hydrogens are added based on definitions in 
the [dic\_lipids.py](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/dic_lipids.py) dictionary. 
This is done by giving a dictionary where the key is the molecule name (abbreviation listed in table above) 
and the value is the correct dictionary key in [dic\_lipids.py](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/dic_lipids.py). 
If correct dictionary key is not yet found, you need to add to [dic\_lipids.py](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/dic_lipids.py). 
In the case of an all atom simulation, UNITEDATOM is left empty.

#### TYPEOFSYSTEM (compulsory)
Tell whether the system is lipid bilayer or something else. Only lipid bilayers are currently supported, but other systems will be included in the future.

#### PUBLICATION
Give reference to a publication(s) related to the data.

#### AUTHORS\_CONTACT
Give the name and email of the main author(s) of the data.

#### SYSTEM
Give description of system in free format. For example ''POPC with cholesterol at 301K''.

#### SOFTWARE\_VERSION
Give the version of the software used.

#### FF
Give the name of the force field used used in the simulation.

#### FF\_SOURCE
Describe the source of the force field parameters. For example, CHARMM-GUI, link to webpage where parameters were downloaded, or citation to a paper.

#### FF\_DATE
Give the date when parameters were accessed or created. The format is day/month/year.

#### Individual force field names for molecules
In some cases special force fields are used for certain molecules. For example, non-standard parameters for ions or other molecules have been used. These can be specified giving forcefield names separately for individual molecules. These can be given as parameters named as
FF+{abbreviation from table above}, i.e., FFPOPC, FFPOT, FFSOL etc.

#### CPT (Gromacs)
Give the name of the Gromacs checkpoint file that is found from the DOI given above.
CPT stands for the name of the cpt file. 

#### LOG (Gromacs)
Give the name of the Gromacs log file that is found from the DOI given above.

#### TOP (Gromacs)
Give the name of the Gromacs top file that is found from the DOI given above.

#### TRAJECTORY\_SIZE
Size of the trajectory file in bytes.

#### TRJLENGTH
Lenght of the trajectory (ps).

#### TEMPERATURE
Temperature of the simulation.

#### NUMBER\_OF\_ATOMS
Total number of atoms in the simulation.

#### DATEOFRUNNIG
Date when added into the databank.

#### EXPERIMENT
Potentially connected experimental data.

#### COMPOSITION (output)
When adding a simulation into NMRlipids databank with the ``AddData.py``, numbers of molecules are automatically calculated and stored into ``COUNT`` keys for each molecule in the COMPOSITION dictionary. Number of lipids are calculated separately for both membrane leaflets. For example, the result for the simulation with COMPOSITION input exemplified above is:


        COMPOSITION:
         CHOL:
          NAME: CHL1
          COUNT:
           - 25
           - 25
          MAPPING: mappingCHOLESTEROLcharmm.yaml
         POPC:
          NAME: POPC
          COUNT:
           - 100
           - 100
          MAPPING: mappingPOPCcharmm.yaml
         SOL:
          NAME: TIP3
          COUNT: 9000
          MAPPING: mappingTIP3PCHARMMgui.yaml
         SOD:
          NAME: SOD
          COUNT: 21
          MAPPING: mappingSOD.yaml
         CLA:
          NAME: CLA
          COUNT: 21
          MAPPING: mappingCLA.yaml
