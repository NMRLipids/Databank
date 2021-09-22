# User input
The databank entries are given as yaml files which locate in subfolders of [Scripts/BuildDatabank/info_files/](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/info_files/) folder with integer names.
The dictionary variables that are used for indexing the data entries are described here.
The necessary parameters required for the analyses are described first and are marked as ''compulsory''. 
Because also other than compulsory parameters may be highly useful in upcycling the data, we strong recommend that all the possible parameters are given upon entry.

When making the entry file for your simulation, looking already existing examples in 
subfolders of [Scripts/BuildDatabank/info_files/](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/info_files/) is useful.

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

#### Molecule names (compulsory)
In the databank, each molecule has a unique abbreviation (see table below). 
You need to give the residue names (in Gromacs: the fourth column of the \[atoms\] directive in an itp-file) of these molecules corresponding your simulation.
If atoms in a lipid belong to different residues (typical situation in Amber force fields), 
give the name of the head group residue here, and add the residue name of each atom to the third column in the mapping file (see below). 
If your simulation contains molecules that are not yet in the databank, you need to define the abbreviation and add molecules 
to the lipids\_dict, molecules\_dict, molecule\_numbers\_dict and molecule\_ff\_dict in 
the [databankLibrary.py](https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/databankLibrary.py), as well as to table below. 

Abbreviation | Molecule name 
------------ | -------------
POPC |  1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine
POPG |  1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoglycerol
POPS | 1-palmitoyl-2-oleoyl-sn-glycero-3-phospho-L-serine
POPE | 1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine
DMPC | 1,2-dimyristoyl-sn-glycero-3-phosphocholine
DPPC | 1,2-dipalmitoyl-sn-glycero-3-phosphocholine
DPPE | 1,2-dipalmitoyl-sn-glycero-3-phosphoethanolamine
DPPG | 1,2-dipalmitoyl-sn-glycero-3-phospho-(1'-rac-glycerol) (sodium salt)
DEPC | 1,2-dierucoyl-sn-glycero-3-phosphocholine
DLPC | 1,2-dilauroyl-sn-glycero-3-phosphocholine
DLIPC| 1,2-dilinoleoyl-sn-glycero-3-phosphocholine
DOPC | 1,2-dioleoyl-sn-glycero-3-phosphocholine
DDOPC| 1,2-didocosahexaenoyl-sn-glycero-3-phosphocholine
DOPS | 1,2-dioleoyl-sn-glycero-3-phospho-L-serine
DSPC | 1,2-distearoyl-sn-glycero-3-phosphocholine
DAPC | 1,2-diarachidonoyl-sn-glycero-3-phosphocholine
POPI | 
SAPI | 
SLPI | 
CHOL | cholesterol 
DHMDMAB | dihexadecyldimethylammonium 
POT | potassium ion 
SOD | sodium ion 
CLA | chloride ion
CAL | calcium ion 
SOL | water 
    
#### MAPPING\_DICT (compulsory)
For the analysis, we need to know the names of atoms in your system.
These are defined using [the mapping file convention](\url{http://nmrlipids.blogspot.com/2015/03/mapping-scheme-for-lipid-atom-names-for.html}), 
where first column gives the universal atom name, second gives the atom name in force field, and third gives the residue name if not the same for all atoms.
The name of the mapping file for each molecule needs to be given in dictionary format. The dictionary the key is 
the molecule name (abbreviation listed in table above) and the value is the name of the mapping file. 
The already existing mapping files can be found from the 
[Scripts/BuildDatabank/mapping_files](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank/mapping_files)
directory. 
If a mapping file for your molecule(s) do not exist, you need to construct one and add to 
the [mapping\_files](https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank/mapping_files) directory.
The mapping file should contain all the atoms of the molecules and no extra lines. 
The number of lines is used to calculate the number of atoms according to the databank which is compared to the number of atoms in trajectory for consistency check.

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
