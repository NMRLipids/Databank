import os
import urllib.request
import yaml

"""
This script is made library-independent by intention. 
We want to be able to run it without additional dependecies.
TODO: remove additional dependencies and import only directroy names
"""

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
    
    def __init__(self, path=r"../../Data/Simulations/"):
        self.path = path
        self.systems = []
        self.__load_systems__(path)
        print('Databank initialized from the folder:', os.path.realpath(path))

    def __load_systems__(self, path):
        for subdir, dirs, files in os.walk(path):
            for filename in files:
                filepath = os.path.join(subdir, filename)
                #print(filepath)
                if filename == "README.yaml":
                    with open(filepath) as yaml_file:
                        content = yaml.load(yaml_file, Loader=yaml.FullLoader)
                        size = len(filepath)
                        sizePath = len(path)
                        content["path"] = filepath[sizePath : size - 11]
                        self.systems.append(content)

    def get_systems(self):
        """ Returns a list of all systems in the NMRlipids databank """
        return self.systems


def download_link(doi, file):  # deprecated?
    """    :meta private:"""
    if "zenodo" in doi.lower():
        zenodo_entry_number = doi.split(".")[2]
        return "https://zenodo.org/record/" + zenodo_entry_number + "/files/" + file
    else:
        print("DOI provided: {0}".format(doi))
        print(
            "Repository not validated. Please upload the data for example to zenodo.org"        )
        return ""

path = '../../Data/Simulations/'
db_data = databank(path)
systems = db_data.get_systems()


for system in systems:
    path = '../../Data/Simulations/' + system['path']
    outfilename = path + 'conf.pdb'
    if os.path.isfile(outfilename):
        continue

    print('Analyzing: ', system['path'])

    try:
        doi = system.get('DOI')
        tpr = system.get('TPR')
        trj_name = path + system.get('TRJ')[0][0]
        tpr_name = path + system.get('TPR')[0][0]
        tpr_url = download_link(doi, tpr[0][0])
        
        if (not os.path.isfile(tpr_name)):
            response = urllib.request.urlretrieve(tpr_url, tpr_name)
    
        if 'gromacs' in system['SOFTWARE']:
            os.system('gmx editconf -f '+ tpr_name + ' -o ' + outfilename)
            
    except:
        pass


        

