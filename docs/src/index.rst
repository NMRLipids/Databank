.. NMRlipids databank documentation master file, created by
   sphinx-quickstart on Mon Sep  4 12:07:48 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _index:

NMRlipids Databank Project
==========================

NMRlipids databank is a community-driven catalogue of biologically relevant lipid membranes emerging from the `NMRlipids Open Collaboration <http://nmrlipids.blogspot.com/>`_. It gathers atomistic MD simulations, X-ray scattering experiments, and NMR experiments on C-H bond order parameters in one place.


*NMRlipids databank is an overlay databank.* Each databank entry (molecule, simulation or experiment) contains the metadata YAML-file, which stores all the essential information for the data upcycling and reuse. This includes, for example, the information about permanent location of each simulation file, but raw data is located in distributed locations outside the NMRlipids Databank. The organisation of the data as well as description of metadata schemas are described in :ref:`dbstructure`. Physically, the database is located in `BilayerData GitHub repository <https://github.com/NMRlipids/BilayerData>`_.

The scientific background and initial motivation of the project is described in the `NMRlipids Databank Manuscript (Nat.Comm., 2024) <https://doi.org/10.1038/s41467-024-45189-z>`_.

databank.nmrlipids.fi web-UI
----------------------------

`NMRlipids Databank-webUI <https://databank.nmrlipids.fi/>`_ provides easy access to the NMRlipids Databank content. Simulations can be searched based on their molecular composition, force field,
temperature, membrane properties, and quality; the search results are ranked based on the simulation quality as evaluated
against experimental data when available. Web-UI provides basic graphical reports
for the computed properties as well as graphical comparison between simulation
and experimental data.

The Web-UI is being developed in the repository `BilayerGUI_laravel <https://github.com/NMRlipids/BilayerGUI_laravel>`_.

NMRlipids Databank-API
----------------------
The NMRlipids Databank-API provides programmatic access to all simulation data in the NMRlipids Databank.
This enables wide range of novel data-driven applications from construction of machine learning models that predict membrane properties, to automatic analysis of virtually
any property across all simulations in the Databank. For examples of novel analyses enabled by the NMRlipids databank API see the `NMRlipids databank manuscript <https://doi.org/10.26434/chemrxiv-2023-jrpwm>`_.

Functions available for simulation analyses are described in :ref:`APIfunctions`.
A project `template <https://github.com/NMRLipids/databank-template>`_ designed to intialize projects that analyse data from NMRlipids databank contains
a `minimum example for looping over available simulations <https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb>`_.
For further examples, see codes that analyze the `area per lipid <https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/calcAPL.py>`_,
`C-H bond order parameters <https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/calcOrderParameters.py>`_,
`X-ray scattering form factors <https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/calc_FormFactors.py>`_,
and `principal component equilibration <https://github.com/NMRLipids/Databank/blob/main/Scripts/AnalyzeDatabank/NMRPCA_timerelax.py>`_.
For these analyses, the universal molecule and atom names are connected to simulation specific names using README.yaml and mapping files
as described in :ref:`molecule_names`.


Adding simulations into the NMRlipids databank
------------------------------------

The NMRlipids Databank is open for additions of simulation and experimental data by anyone. For detailed instructions to add new data, to update databank analyses and run quality evaluations, see :ref:`addSimulation`. Quick and minimal steps to add a new simulation are here:

#. Add trajectory and topology (tpr for Gromacs, pdb or corresponding to other programs) file into a `Zenodo <https://zenodo.org/>`_ repository.

#. Login to `NMRlipids Upload Portal <https://upload-portal.nmrlipids.fi/>`_ and fill up the form with metadata for your simulation. For instructions, see :ref:`readmesimu` and `examples <https://github.com/NMRLipids/BilayerData/tree/main/info_files>`_. Mapping files are described in  :ref:`molecule_names` and are available from `here <https://github.com/NMRLipids/BilayerData/tree/main/Molecules/membranes>`_ .

#. Your simulation will be automatically processed via GitHub Actions on the server side after the approval of your submission by one of the NMRlipids contributors.

Experimental data addition is currently not automatized and should be performed manually via making pull request to the BilayerData repository. The instrutions are available at :ref:`addingExpData`.

Do not hesitate to ask assistance regarding data addition on the GitHub page of `BilayerData Issues <https://github.com/NMRLipids/BilayerData/issues>`_.


Installation and system requirements
------------------------------------

The code has been tested in Linux, MacOS and Windows environment with python 3.10 or newer and recent `Gromacs <https://manual.gromacs.org/current/install-guide/index.html>`_ version installed.

Setup is straingforward using ``pip`` or ``uv pip``:

.. code-block:: bash

   pip install nmrlipids_databank
   nml_initialize_data dev
   source databank_env.rc

More detailed instructions are coming soon.

**TODO: write the block**



.. toctree::
   :maxdepth: 3
   :caption: Python Interface

   gettingstarted
   dbprograms
   Overview

.. toctree::
   :maxdepth: 1
   :caption: Membrane Databank

   dbstructure
   dbcontribute

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
