import os, sys
import urllib.request

"""
Creates PDB for every GROMACS-system and 
It imports just `core` and `databankio` to avoid additional dependecies.
"""
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from DatabankLib.core import initialize_databank
from DatabankLib import NMLDB_SIMU_PATH
from DatabankLib.databankio import resolve_download_file_url

systems = initialize_databank()

for system in systems:
    path =  os.path.join(NMLDB_SIMU_PATH, system['path'])
    outfilename = os.path.join(path, 'conf.pdb')
    if os.path.isfile(outfilename):
        continue

    print('Analyzing: ', system['path'])

    try:
        doi = system.get('DOI')
        tpr = system.get('TPR')
        tpr_name = os.path.join(path, system.get('TPR')[0][0])
        
        if (not os.path.isfile(tpr_name)):
            tpr_url = resolve_download_file_url(doi, tpr[0][0])
            response = urllib.request.urlretrieve(tpr_url, tpr_name)
    
        if 'gromacs' in system['SOFTWARE']:
            os.system('gmx editconf -f '+ tpr_name + ' -o ' + outfilename)
            
    except:
        pass

