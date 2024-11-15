"""
Core databank class and system initialization function.
Imported by `databankLibrary` by default.
Can be imported without additional libraries to scan Databank system file tree!
"""

import os
import yaml
import collections.abc
from DatabankLib import NMLDB_SIMU_PATH


class SystemsCollection(collections.abc.Sequence):
    """Immutable collection of system dicts. Can be accessed by ID using loc()."""

    def __init__(self, iterable=[]):
        self.data = iterable
        self.__genIndexByID()

    def __genIndexByID(self):
        self._idx = dict()
        for i in range(len(self)):
            if 'ID' in self[i].keys():
                self._idx[self[i]['ID']] = i

    def __getitem__(self, i):
        return self.data[i]

    def __len__(self):
        return len(self.data)

    def loc(self, id: int):
        return self.data[self._idx[id]]


class databank:
    """ :meta private:
    Representation of all simulation in the NMR lipids databank.

        `path` should be the local location of /Data/Simulations/ in the NMRlipids
        databank folder. Example usage to loop over systems:

            path = '../../Data/Simulations/'
            db_data = databank(path)
            systems = db_data.get_systems()

            for system in systems:
                print(system)
    """

    def __init__(self):
        self.path = NMLDB_SIMU_PATH
        __systems = self.__load_systems__()
        self._systems = SystemsCollection(__systems)
        print('Databank initialized from the folder:', os.path.realpath(self.path))

    def __load_systems__(self):
        systems = []
        rpath = os.path.realpath(self.path)
        for subdir, dirs, files in os.walk(rpath):
            for filename in files:
                filepath = os.path.join(subdir, filename)
                if filename == "README.yaml":
                    with open(filepath) as yaml_file:
                        content = yaml.load(yaml_file, Loader=yaml.FullLoader)
                        relpath = os.path.relpath(filepath, rpath)
                        content["path"] = relpath[:-11]
                        systems.append(content)
        return systems

    def get_systems(self):
        """ Returns a list of all systems in the NMRlipids databank """
        return self._systems


def initialize_databank():
    """
    Intializes the NMRlipids databank.

    :return: list of dictionaries that contain the content of README.yaml files for
             each system.
    """
    db_data = databank()
    return db_data.get_systems()


# TODO: is not used at all in the project!!
def print_README(system):
    """
    Prints the content of ``system`` dictionary in human readable format.

    :param system: NMRlipids databank dictionary defining a simulation.

    """
    if system == 'example':
        readmePath = '../data/simulations/READMEexplanations.yaml'  # TODO: CORRECT!
        with open(readmePath, 'r') as file:
            readmeFile = yaml.safe_load(file)
    else:
        readmeFile = system

    for key in readmeFile:
        print('\033[1m' + key + ":" + '\033[0m')
        print(" ", readmeFile[key])
