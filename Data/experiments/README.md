Experimental C-H bond order parameters from NMR (OrderParameters) and x-ray scattering form factors (FormFactors) 
are collected to folder named according to the DOI of the original publication. Experimental conditions are described in dictionary
format in README.yaml files within these folders. The keys of the dictionaries are defined as:

|key | description |
| ----| -----|
|DOI | DOI of the of the original publication of the experimental data.|
TEMPERATURE | Temperature of the experiment.
MOLAR_FRACTIONS | Dictionary of molar fractions of bilayer components.
ION_CONCENTRATIONS | Dictionary of ion concentrations of the system ([ion]= N<sub>ion</sub>*[water]/N<sub>water</sub>, [water] = 55.5 M)
TOTAL_LIPID_CONCENTRATION | Total concentration of lipid components. If exact concentration is not known, but experiments are performed in excess water, ’full hydration’ can be given.
COUNTER_IONS | Type of counter ions if present.
