.. _listOfFiles:

List and descriptions of NMRlipids databank files
=======================================================

README.yaml
-----------
A python dictionary containing all the relevant information of each simulation entry. These files are stored at `/Data/Simulations <https://github.com/NMRLipids/Databank/tree/main/Data/Simulations>`_. For more information see :ref:`readmecontent`.

info.yaml
--------
Contains information given by contributor when adding data into the NMRlipids databank. This file is given as an input to ``AddData.py`` to create ``README.yaml`` files: :code:`python3 AddData.py -f info.yaml`. These files are stored at `BuildDatabank/info_files <https://github.com/NMRLipids/Databank/tree/main/Scripts/BuildDatabank/info_files>`_. For more information see :ref:`readmecontent` and :ref:`addData`.  

mapping files
------------
Describe connection between universal atom names and simulation specific atom names. These enable automatic analyses over simulations with different naming conventions. For more information, see :ref:`molecule_names`

apl.json
--------
Area per lipid (Å\ :sup:`2`) as a function of time (ps) for a simulation. Calculated with ``calcAPL.py`` and stored in `/Data/Simulations <https://github.com/NMRLipids/Databank/tree/main/Data/Simulations>`_ for all simulations.

{lipid name}_OrderParameters.json
---------
C-H bond order parameters calculated with ``calcOrderParameters.py`` and stored in `/Data/Simulations <https://github.com/NMRLipids/Databank/tree/main/Data/Simulations>`_ for all simulations. Key in json format gives the universal C and H atom names. Values are average over lipids, standard deviation, and the standard error of the mean, respectively.

FormFactor.json
---------------
X-ray scattering calculated with ``calc_FormFactors.py`` and stored in `/Data/Simulations <https://github.com/NMRLipids/Databank/tree/main/Data/Simulations>`_ for all simulations. X-axis value unit are Å\ :sup:`-1` and y-axis value units are e/nm\ :sup:`2`.

TotalDensity.json
-----------------
Total electron densities calculated with ``calc_FormFactors.py`` and stored in `/Data/Simulations <https://github.com/NMRLipids/Databank/tree/main/Data/Simulations>`_ for all simulations. X-axis units are nm and y-axis units e/nm\ :sup:`3`.

thickness.json
-------------
Membrane thickness (nm) calculated with ``calc_thickness.py`` and stored in `/Data/Simulations <https://github.com/NMRLipids/Databank/tree/main/Data/Simulations>`_ for all simulations.

eq_times.json
-------------
Relative equilibration time of molecule principal components with ``NMRPCA_timerelax.py`` and stored in `/Data/Simulations <https://github.com/NMRLipids/Databank/tree/main/Data/Simulations>`_ for all simulations. 

{lipid name}_OrderParameters_quality.json
---------------------------
Dictionary with the quality of each C-H bond against experiments if available. First key is the DOI for the source of experimental data. Second key gives C-H pair univeral atom names. Third key gives values for simulation order parameter, its stardard deviation, standard error of the mean, experimental order parameter, its error, and finally the quality of the order parameter. Quality is the probability for the agreement between simulated and experimental values taking into account the error bars. For more details, see the `NMRlipids databank manuscript <https://doi.org/10.26434/chemrxiv-2023-jrpwm>`_.



{lipid name}_FragmentQuality.json
--------------------------------
Fragment qualities determined separately for each lipid in the simulation with experimental data available using Eq. (4) in the `NMRlipids databank manuscript <https://doi.org/10.26434/chemrxiv-2023-jrpwm>`_.


SYSTEM_quality.json
-------------------
Total qualities averaged over different lipids for different membrane parts calculated from Eq. (5) in the `NMRlipids databank manuscript <https://doi.org/10.26434/chemrxiv-2023-jrpwm>`_.

FormFactorQuality.json
----------------------
Quality of form factor against experimental data determined as described in the `NMRlipids databank manuscript <https://doi.org/10.26434/chemrxiv-2023-jrpwm>`_. Second term is the scaling coefficient for experimental intensities (Eq. (6) in the `NMRlipids databank manuscript <https://doi.org/10.26434/chemrxiv-2023-jrpwm>`_).

