"""
Core databank class and system initialization function.
Imported by `databankLibrary` by default.
Can be imported without additional libraries to scan Databank system file tree!
"""

import os
import yaml
from . import NMLDB_SIMU_PATH

class databank:
    """ :meta private: 
    Representation of all simulation in the NMR lipids databank. 

        `path` should be the local location of /Data/Simulations/ in the NMRlipids databank folder. Example usage to loop over systems: 
   
            path = '../../Data/Simulations/'
            db_data = databank(path)
            systems = db_data.get_systems()

            for system in systems:
                print(system)
    """
    
    def __init__(self):
        self.path = NMLDB_SIMU_PATH
        self.systems = []
        self.__load_systems__()
        print('Databank initialized from the folder:', os.path.realpath(self.path))

    def __load_systems__(self):
        rpath = os.path.realpath(self.path)
        for subdir, dirs, files in os.walk(rpath):
            for filename in files:
                filepath = os.path.join(subdir, filename)
                if filename == "README.yaml":
                    with open(filepath) as yaml_file:
                        content = yaml.load(yaml_file, Loader=yaml.FullLoader)
                        relpath = os.path.relpath(filepath, rpath)
                        content["path"] = relpath[ : -11]
                        self.systems.append(content)

    def get_systems(self):
        """ Returns a list of all systems in the NMRlipids databank """
        return self.systems



def initialize_databank():
    """ 
    Intializes the NMRlipids databank.

    :return: list of dictionaries that contain the content of README.yaml files for each system.  
    """
    db_data = databank()
    systems = db_data.get_systems()
    return systems

#TODO: is not used at all in the project!!
def print_README(system):
    """ 
    Prints the content of ``system`` dictionary in human readable format. 

    :param system: NMRlipids databank dictionary defining a simulation.

    """
    if system == 'example':
        readmePath = '../data/simulations/READMEexplanations.yaml' #TODO: CORRECT!
        with open(readmePath, 'r') as file:
            readmeFile = yaml.safe_load(file)
    else:
        readmeFile = system
            
    for key in readmeFile:
        print('\033[1m' + key + ":" + '\033[0m')
        print(" ", readmeFile[key])

