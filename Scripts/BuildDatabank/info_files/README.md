# User input and content of README.yaml files
README.yaml files in subfolders of [Data/Simulations/](https://github.com/NMRLipids/Databank/tree/main/Data/Simulations) contain all the essential information of the simulations in the NMRlipids databank. 
These files are build based on information given in files locating in subfolders of [Scripts/BuildDatabank/info_files/](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/info_files/). 
Each entry is in subfolder with integer name.
The manually (compulsory and optional) entered and automatically exctracted information in the README.yaml files are listed in the table below.
The necessary parameters required for the analyses are compulsory, but addition of also optional parameters is highly recommended whenever possible.

The already existing info files in 
subfolders of [Scripts/BuildDatabank/info_files/](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/info_files/)
serve as useful examples when preparing the new entry files.

key | description | type  
----|--------|------
DOI | DOI from where the raw data is found | user given (compulsory) 
TRJ | Name of the trajectory file found from DOI | 
TPR | Name of the topology file found from DOI (tpr file in the case of Gromacs) | 
SOFTWARE | Software used to run the simulation (e.g. Gromacs, Amber, NAMD, etc.) | 
PREEQTIME | Pre-equilibrate time simulated before the uploaded trajectory in nanoseconds. For example, if you upload 100-200 ns part of total 200 ns simulation, this should value should be 100. | 
TIMELEFTOUT | Equilibration period in the uploaded trajectory that should be discarded in analyses. For example, if you upload 0-200 ns part of total 200 ns simulation where the first 100 ns should be considered as an equilibration, this value should be 100. 
COMPOSITION | Molecules names used in the simulation and corresponding mapping files. For more detailed description see below. | 
DIR\_WRK | Temporary working directory in your local computer. 
UNITEDATOM\_DICT | Information for constructing hydrogens for united atom simulations, empty for all atom simulations | 
TYPEOFSYSTEM | Lipid bilayer or something else | 
||
PUBLICATION | Give reference to a publication(s) related to the data. | User given (optional)
AUTHORS\_CONTACT | Name and email of the main author(s) of the data. | 
SYSTEM | System description in the free text format | 
SOFTWARE\_VERSION | Version of the used software | 
FF | Name of the used force field | 
FF\_SOURCE | Source of the force field parameters, e.g, CHARMM-GUI, webpage, citation to a publication, etc. | 
FF\_DATE |  Date when force field parameters were accessed in the given source (day/month/year). | 
FF{molename} | Molecule specific force field information, e.g., water model with FFSOL and sodium parameters with FFSOD. | 
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
COMPOSITION | Numbers of lipid molecules (NPOPC, NPOPG, etc.) per membrane leaflet are calculated by determining on which side of the center of mass of the membrane the center of mass of the head group of each lipid molecule is located. Numbers of other molecules such as solvent and ions (NSOL, NPOT, NSOD, etc.) are read from the topology file. | 

#### DOI (compulsory)
Give the DOI identity for the location of simulation files. 
Current databank works only for the data in [Zenodo](www.zenodo.org), but other potential sources are may be implemented in the future. 
Note that the DOI must point to a specific version of dataset in Zenodo. DOIs pointing to all versions of certain dataset do not work.

#### SOFTWARE (compulsory)
Give the name of software used to run the simulation. The options are GROMACS, AMBER, NAMD, CHARMM and OPENMM. So far, only simulations run with GROMACS are accepted by the script.

#### TRJ (compulsory)
Give the name of the trajectory file that is found from the DOI given above.

#### TPR (compulsory)
Give the name of the file with topology information (tpr file in the case of Gromacs) that is found from the DOI given above.

#### PREEQTIME (compulsory)
Give the time simulated before the uploaded trajectory in nanoseconds. For example, if you upload 100-200 ns part of total 200 ns simulation, this should value should be 100.

#### TIMELEFTOUT (compulsory)
Give the time that should be considered as an equilibration period in the uploaded trajectory. Frames before the give time will be discarded in the analysis. 
For example, if you upload 0-200 ns part of total 200 ns simulation where the first 100 ns should be considered as an equilibration, this value should be 100.

#### COMPOSITION (compulsory)
Information that is used to determine the composition (i.e., the number of molecules) of each simulation is given in python dictionary format in COMPOSITION.
The COMPOSITION is given in a dictionary format, where the first key is the universal molecule name (abbreviation listed in the table below).
Another dictionary is then created for each molecule containing the values for the residue name in your simulation (NAME) and the name of the mapping file (MAPPING). 
For example, the COMPOSITION dictionary for [a system](https://doi.org/10.5281/zenodo.259392) containing POPC, Cholesterol (CHOL), water (SOL), sodium (SOD), and chloride (CLA) can be given as

COMPOSITION:<br />
&nbsp;&nbsp;CHOL:<br />
&nbsp;&nbsp;&nbsp;&nbsp;NAME: CHL1<br />
&nbsp;&nbsp;&nbsp;&nbsp;MAPPING: mappingCHOLESTEROLcharmm.txt<br />
&nbsp;&nbsp;POPC:<br />
&nbsp;&nbsp;&nbsp;&nbsp;NAME: POPC<br />
&nbsp;&nbsp;&nbsp;&nbsp;MAPPING: mappingPOPCcharmm.txt<br />
&nbsp;&nbsp;SOL:<br />
&nbsp;&nbsp;&nbsp;&nbsp;NAME: TIP3<br />
&nbsp;&nbsp;&nbsp;&nbsp;MAPPING: mappingTIP3PCHARMMgui.txt<br />
&nbsp;&nbsp;SOD:<br />
&nbsp;&nbsp;&nbsp;&nbsp;NAME: SOD<br />
&nbsp;&nbsp;&nbsp;&nbsp;MAPPING: mappingSOD.txt<br />
&nbsp;&nbsp;CLA:<br />
&nbsp;&nbsp;&nbsp;&nbsp;NAME: CLA<br />
&nbsp;&nbsp;&nbsp;&nbsp;MAPPING: mappingCLA.txt<br />

For example in Gromacs, the residue name can be found from the fourth column of the \[atoms\] directive in an itp-file.
Mapping files define the atom names in your system using [the mapping file convention](\url{http://nmrlipids.blogspot.com/2015/03/mapping-scheme-for-lipid-atom-names-for.html}), 
where first column gives the universal atom name, second gives the atom name in force field, and third gives the residue name if not the same for all atoms.
The already existing mapping files can be found from the 
[Scripts/BuildDatabank/mapping_files](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank/mapping_files)
directory. 
If a mapping file for your molecule(s) do not exist, you need to construct one and add to 
the [mapping\_files](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank/mapping_files) directory.
If atoms in a lipid belong to different residues (typical situation in Amber force fields), 
give the name of the head group residue in COMPOSITION dictionary, and add the residue name of each atom to the third column in the mapping file. 
If your simulation contains molecules that are not yet in the databank, you need to define the abbreviation and add molecules 
to the lipids\_dict, molecules\_dict, molecule\_numbers\_dict and molecule\_ff\_dict in 
the [databankLibrary.py](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/databankLibrary.py), as well as to table below. 
The mapping file should contain all the atoms of the molecules and no extra lines. 
The number of lines is used to calculate the number of atoms according to the databank which is compared to the number of atoms in trajectory for consistency check.

Abbreviation | Molecule name 
------------ | -------------
POPC |  1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine
POPG |  1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoglycerol
POPS | 1-palmitoyl-2-oleoyl-sn-glycero-3-phospho-L-serine
POPE | 1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine
PYPC | 1-(16:0)-2-(16:1$^\Delta9$)-sn-glycero-3-phosphocholine
PAzePCprot | 1-palmitoyl-2-azelaoyl-sn-glycero-3-phosphocholine protonated
PAzePCdeprot | 1-palmitoyl-2-azelaoyl-sn-glycero-3-phosphocholine deprotonated
DMPC | 1,2-dimyristoyl-sn-glycero-3-phosphocholine
DPPC | 1,2-dipalmitoyl-sn-glycero-3-phosphocholine
DPPE | 1,2-dipalmitoyl-sn-glycero-3-phosphoethanolamine
DPPG | 1,2-dipalmitoyl-sn-glycero-3-phospho-(1'-rac-glycerol) (sodium salt)
DEPC | 1,2-dierucoyl-sn-glycero-3-phosphocholine
DRPC | 1,2-(14:1$^\Delta9$)-sn-glycero-3-phosphocholine
DYPC | 1,2-(16:1$^\Delta9$)-sn-glycero-3-phosphocholine
DLPC | 1,2-dilauroyl-sn-glycero-3-phosphocholine
DLIPC| 1,2-dilinoleoyl-sn-glycero-3-phosphocholine
DOG  | 1,2-dioleoyl-sn-glycerol
DOPC | 1,2-dioleoyl-sn-glycero-3-phosphocholine
DOPE | 1,2-dioleoyl-sn-glycero-3-phosphoethanolamine
DDOPC| 1,2-didocosahexaenoyl-sn-glycero-3-phosphocholine
DOPS | 1,2-dioleoyl-sn-glycero-3-phospho-L-serine
DSPC | 1,2-distearoyl-sn-glycero-3-phosphocholine
DAPC | 1,2-diarachidonoyl-sn-glycero-3-phosphocholine
SLiPC | 1-(18:0)-2-(18:2 $^{\Delta9,12}$)-sn-glycero-3-phosphocholine 
DMTAP | 1,2-dimyristoyl-3-trimethylammonium-propane
SOPC | 1-stearoyl-2-oleoyl-sn-glycero-3-phosphocholine
POPI | 
SAPI | 
SLPI | 
SDG | 1-stearoyl-2-docosahexaenoyl-sn-glycerol
SDPE | 1-stearoyl-2-docosahexaenoyl-sn-glycero-3-phosphoethanolamine
CER  | N-palmitoyl-D-erythro-sphingosine
CHOL | cholesterol 
DCHOL | 18,19-di-nor-cholesterol
DHMDMAB | dihexadecyldimethylammonium 
POT | potassium ion 
SOD | sodium ion 
CLA | chloride ion
CAL | calcium ion 
CES | caesium ion
C20 | n-eicosane
SOL | water 
    
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
