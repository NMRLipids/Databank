.. _APIfunctions:

NMRlipids databank API functions
======================

The functions that help analyzing data in the NMRlipids databank are described here. These functions are located in `Databank/Scripts/DatabankLib/databankLibrary.py <https://github.com/NMRLipids/Databank/blob/main/Scripts/DatabankLib/databankLibrary.py>`_. To get started using these functions, first create the folder where you want to work

.. code-block::

   mkdir NMRlipids/
   cd NMRlipids
   
then clone the NMRlipids databank git into this folder

.. code-block::

   git clone https://github.com/NMRLipids/Databank.git
   cd Databank
   pip install -e .

You can then start to work with the `template <https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb>`_ or write a code from the scratch. The minimum python code to intialize NMRlipids databank is

.. code-block::

   import DatabankLib
   from DatabankLib.core import initialize_databank

   systems = initialize_databank()

After running this, ``systems`` is a SystemsCollection which works like a list but with added functionality and contains dictionaries where each dictionary is a simulation in the NMRlipids databank. A simulation dictionary contains the content of the README.yaml for that simulation. The content of README.yaml files is described in :ref:`readmecontent`. ``systems`` can be then used to loop over all simulations:

.. code-block::

   for system in systems:
       print(system)
   
Examples on analyses over NMRlipids databank can be found from the `template <https://github.com/NMRLipids/databank-template/blob/main/scripts/template.ipynb>`_ and `codes used to analyze the results for the NMRlipids databank manuscript <https://github.com/NMRLipids/DataBankManuscript/tree/main/scripts>`_.


The functions available to analyze the simulations can be found :doc:`here <auto_gen/Scripts.DatabankLib.databankLibrary>`

More information on initializing the databank and related information can be found on the  :doc:`DatabankLib.core page. <auto_gen/Scripts.DatabankLib.core>` 

