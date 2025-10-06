.. _gettingstarted:

Getting started
===============

The functions that help analyzing data in the NMRlipids databank are described here. These functions are located in the :mod:`DatabankLib` module. To get started using these functions, first set up the package and initialize the databank:

.. code-block:: bash

   pip install nmrlipids_databank
   nml_initialize_data toy
   source databank_env.rc

"Toy" is a small test databank which is useful for testing and learning the package. You can then start to work with the `template <https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb>`_ or write a code from the scratch. The minimum python code to intialize NMRlipids databank is

.. code-block:: python

   import DatabankLib
   from DatabankLib.core import initialize_databank

   systems = initialize_databank()

After running this, ``systems`` is an instance of :class:`DatabankLib.core.SystemsCollection` which works like a list but with added functionality and contains dictionaries where each dictionary is a simulation in the NMRlipids databank.
A simulation dictionary contains the content of the README.yaml for that simulation.
The content of README.yaml files is described in :ref:`Simulation metadata <readmesimu>`. ``systems`` can be then used to loop over all simulations:

.. code-block:: python

   for system in systems:
       print(system)

Examples on analyses over NMRlipids databank can be found from the `template <https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb>`_ and `codes used to analyze the results for the NMRlipids databank manuscript <https://github.com/NMRLipids/DataBankManuscript/tree/main/scripts>`_.


The core functionality of the DatabankLib including initializing the databank and the definition of main data models is provided by :mod:`DatabankLib.core`. Some useful functions for analysis of the simulations can be found in the :mod:`DatabankLib.databankLibrary` module.


.. toctree::
   :maxdepth: 1

   api/exampleAndTutorials
