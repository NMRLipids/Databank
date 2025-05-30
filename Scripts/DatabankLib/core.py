"""
Core databank class and system initialization function.
Imported by `databankLibrary` by default.
Can be imported without additional libraries to scan Databank system file tree!
"""

import os
import yaml
import collections.abc
from typing import Dict
from DatabankLib.settings.molecules import Molecule

from DatabankLib import NMLDB_SIMU_PATH
from DatabankLib.settings.molecules import Lipid, lipids_set, lipids_dict, molecules_set, NonLipid


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


class System(collections.abc.MutableMapping):
    """
    Main Databank single object, which is an extension of a
    dictionary with additional functionality.
    """

    def __init__(self, data=None):
        self._store = {}
        if isinstance(data, dict):
            self._store.update(data)
        elif isinstance(data, collections.abc.MutableMapping):
            self._store.update(dict(data))
        else:
            raise TypeError("Expected dict or Mapping")

        self._content = {}
        for k,v in self["COMPOSITION"].items():
            mol = None
            if k in lipids_set:
                mol = Lipid(k)
            elif k in molecules_set:
                mol = NonLipid(k)
            mol.register_mapping(v["MAPPING"])
            self._content[k] = mol

    def __getitem__(self, key):
        return self._store[key]

    def __setitem__(self, key, value):
        self._store[key] = value

    def __delitem__(self, key):
        del self._store[key]

    def __iter__(self):
        return iter(self._store)

    def __len__(self):
        return len(self._store)

    @property
    def content(self) -> Dict[str, Molecule]:
        """ Returns dictionary of molecule objects. """
        return self._content


class Databank:
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
                        content = System(yaml.load(yaml_file, Loader=yaml.FullLoader))
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
    db_data = Databank()
    return db_data.get_systems()


# TODO: is not used at all in the project!!
def print_README(system):
    """
    Prints the content of ``system`` dictionary in human readable format.

    :param system: NMRlipids databank dictionary defining a simulation.

    """
    if system == 'example':
        current_folder = os.path.dirname(os.path.realpath(__file__))
        readmePath = os.path.join(current_folder, 'settings', 'READMEexplanations.yaml')
        with open(readmePath, 'r') as file:
            readmeFile = yaml.safe_load(file)
    else:
        readmeFile = system

    for key in readmeFile:
        print('\033[1m' + key + ":" + '\033[0m')
        print(" ", readmeFile[key])
