(addingExpData)=

Adding experimental data into the NMRlipids databank
==================

Experimental data is stored in [Data/experiments](https://github.com/NMRLipids/Databank/tree/main/Data/experiments) folder. C-H bond order parameters from NMR are in [OrderParameters subfolder](https://github.com/NMRLipids/Databank/tree/main/Data/experiments/OrderParameters) and X-ray scattering form factors in [FormFactors subfolder](https://github.com/NMRLipids/Databank/tree/main/Data/experiments/FormFactors). The data is located in subfolders named after the DOI of the original source publication of the data. Because a publication can contain more than one experimental data set, each data set is stored in a subfolder with integer name, for example, [Data/experiments/OrderParameters/10.1039/c2cp42738a/2](https://github.com/NMRLipids/Databank/tree/main/Data/experiments/OrderParameters/10.1039/c2cp42738a/2). Each of such folders then contain the experimental data and `README.yaml` files describing the experimental conditions. The keys of these dictionaries are summarized in the table and described more detailed below.

|key | description |
| ----| -----|
|DOI | DOI of the of the original publication of the experimental data|
TEMPERATURE | Temperature of the experiment
MOLAR_FRACTIONS | Dictionary of molar fractions of bilayer components
ION_CONCENTRATIONS | Dictionary of ion concentrations in the system
TOTAL_LIPID_CONCENTRATION | Total concentration of lipid components
COUNTER_IONS | Type of counter ions if present

DOI:
----------
DOI of the original publication where the experimental data originates.

TEMPERATURE
------------
Temperature of the experiment.

MOLAR_FRACTIONS
-----------------
Dictionary of molar fractions of bilayer components. For example:
```
MOLAR_FRACTIONS:
  POPC: 0.93
  CHOL: 0.07
  ````


ION_CONCENTRATIONS
--------------

Dictionary of ion concentrations of the system ([ion]= N<sub>ion</sub>*[water]/N<sub>water</sub>, where [water] = 55.5 M). For example:
```
ION_CONCENTRATIONS:
  POT: 0
  SOD: 0.1
  CLA: 0.1
  CAL: 0
```


TOTAL_LIPID_CONCENTRATION
---------------------

Total concentration of lipid components ([lipid]= N<sub>lipid</sub>*[water]/N<sub>water</sub>, where [water] = 55.5 M). If exact concentration is not known, but experiments are performed in excess water, ’full hydration’ can be given.

COUNTER_IONS
-------------

Type of counter ions if present.

Steps to add experimental data
-----------------------------

1. Create the above described `README.yaml` file of your data.

2. Copy this README file data into a appropriate directory named as described above.

3. If you have order parameter data, create a file named `{lipidname}_Order_Parameters.dat` where `{lipidname}` is the universal name of the lipid from which the data is measured from. The first two columns of this file should define the atom pair with universal atom names, third column has the experimental order parameter value, and fourth column has the experimental error. If the experimental error is not known, set it to 0.02. Store the created `{lipidname}_Order_Parameters.dat` file into the appropriate folder with the `README.yaml` file. Create `json` version from `dat` file by running
   ```
   python data_to_json.py {dat file path}
   ```
   in folder [Data/experiments/OrderParameters](https://github.com/NMRLipids/Databank/tree/main/Data/experiments/OrderParameters). For example, see subfolders of [Data/experiments/OrderParameters](https://github.com/NMRLipids/Databank/tree/main/Data/experiments/OrderParameters)

4. If you have X-ray scattering form factor data, store the form factor into appropriate folders in ASCII format where first column in x-axis values (Å<sup>-1</sup>), second column is y-axis value, and third is the error. Then create `json` file by runnig
   ```
   python data_to_json.py {ASCII file path}
   ```
   in folder [Data/experiments/FormFactors](https://github.com/NMRLipids/Databank/tree/main/Data/experiments/FormFactors). For example, see subfolders of [Data/experiments/FormFactors](https://github.com/NMRLipids/Databank/tree/main/Data/experiments/FormFactors)

5. To update the connections with simulations after adding new experimental data, perform the step 6 in [Adding simulations into the NMRlipids databank](addData).