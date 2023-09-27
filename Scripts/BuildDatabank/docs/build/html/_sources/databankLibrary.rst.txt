.. _APIfunctions:

NMRlipids databank API functions
======================

The functions that help analyzing data in the NMRlipids databank are described here. These functions locate in `databanklibrary.py <https://github.com/NMRLipids/Databank/blob/main/Scripts/BuildDatabank/databankLibrary.py>`_. To get started using these functions, first create the folder where you want to work

.. code-block::

   mkdir NMRlipids/
   cd NMRlipids
   
then clone the NMRlipids databank git into this folder

.. code-block::

   git clone https://github.com/NMRLipids/Databank.git

You can then start to work with the `template <https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb>`_ or write a code from the scratch. The minimum python code to intialize NMRlipids databank is

.. code-block::

   import sys
   databankPath =  './Databank/'   # this is the local path for the cloned Databank git
   sys.path.insert(1, databankPath + '/Scripts/BuildDatabank/')
   from databankLibrary import * 
   systems = initialize_databank(databankPath)

After running this, ``systems`` is the list of dictionaries where each dictionary is a simulation in the NMRlipids databank. A simulation dictionary contains the content of the README.yaml for that simulation. The content of README.yaml files is described in :ref:`readmecontent`. ``systems`` can be then used to loop over all simulations:

.. code-block::

   for system in systems:
       print(system)
   
Examples on analyses over NMRlipids databank can be found from the `template <https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb>`_ and `codes used to analyze the results for the NMRlipids databank manuscript <https://github.com/NMRLipids/DataBankManuscript/tree/main/scripts>`_.


The functions available to analyze the simulations are listed in here:

.. automodule:: databankLibrary
   :members:
   :undoc-members:
   :show-inheritance:
