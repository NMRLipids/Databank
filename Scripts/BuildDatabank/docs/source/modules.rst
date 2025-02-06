Project structure
=============

Project has four subfolders:

DatabankLib
-----------

Main project package. Contains most of project's logic and API.

AnalyzeDatabank
---------------

Scripts performing analysis of the entire databank.

BuildDatabank
-------------

.. toctree::
   :maxdepth: 1

   AddData
   QualityEvaluation
   searchDATABANK
   makeRanking
   check_mappings
   buildH_calcOP_test


updateGUI
---------

.. toctree::
   :maxdepth: 1

   createPDBs

tests
-----

Contains pytest-based unit-tests and regression tests over toy-databank.