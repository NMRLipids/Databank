import os
import sys
import urllib.request
from importlib import import_module

"""
Creates PDB for every GROMACS-system and
It imports just `core` and `databankio` to avoid additional dependecies.
It DOES NOT require the package to be pre-installed. Only yaml!
"""
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
dbl = import_module("DatabankLib")
core = import_module("DatabankLib.core")
dbio = import_module("DatabankLib.databankio")
sys.path.pop(0)

if __name__ == "__main__":
    systems = core.initialize_databank()

    for system in systems:
        path = os.path.join(dbl.NMLDB_SIMU_PATH, system['path'])
        outfilename = os.path.join(path, 'conf.pdb')
        if os.path.isfile(outfilename):
            continue

        print('Analyzing: ', system['path'])

        doi = system.get('DOI')

        tpr_name = None
        try:
            tpr = system.get('TPR')
            tpr_name = os.path.join(path, system.get('TPR')[0][0])

            if not os.path.isfile(tpr_name):
                tpr_url = dbio.resolve_download_file_url(doi, tpr[0][0])
                response = urllib.request.urlretrieve(tpr_url, tpr_name)
        except Exception:
            print("There is no TPR for this system. We will try GRO.")
        else:
            os.system('gmx editconf -f ' + tpr_name + ' -o ' + outfilename)
            continue

        gro_name = None
        try:
            gro = system.get('GRO')
            gro_name = os.path.join(path, system.get('GRO')[0][0])

            if not os.path.isfile(gro_name):
                gro_url = dbio.resolve_download_file_url(doi, gro[0][0])
                response = urllib.request.urlretrieve(gro_url, gro_name)
        except Exception:
            print("There is no GRO for this system. We will skip it.")
        else:
            os.system('gmx editconf -f ' + gro_name + ' -o ' + outfilename)
            continue
