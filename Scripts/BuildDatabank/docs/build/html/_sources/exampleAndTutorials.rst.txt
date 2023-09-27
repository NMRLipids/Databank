Examples and tutorials
======================

A template project for analyses utilizing the NMRlipids databank is available in `here <https://github.com/NMRLipids/databank-template/tree/main>`_ where examples and templates for analyses are available from the `scripts <https://github.com/NMRLipids/databank-template/tree/main/scripts>`_ folder. 

`Plotting basic simulation properties <https://github.com/NMRLipids/databank-template/blob/main/scripts/plotSimulation.ipynb>`_
----------------------------------
Plots the basic properties of simulation selected based on its NMRlipids databank ID number.

`Show ranking tables of simulations based in their quality against experimental data <https://github.com/NMRLipids/databank-template/blob/main/scripts/plotQuality.ipynb>`_
---------------
Shows different kinds of rankings of simulations against experimental data.

`Template for more advance API usage <https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb>`_
------------------
Demonstrates the usage of API by three examples. 1) Selects a random simulation and prints the related databank content in human readable format. 2) Shows the readily analyzed properties for the selected random simulation (area per lipid, membrane thickness, relative equilibration times, X-ray scattering form factors, and C-H bond order parameters). 3) Selects a random simulation with the trajectory size below 100Mb and calculates P-N vector angle with respect to membrane normal for all lipids for which P and N atoms are available in headgroup.
