#!/usr/bin/env python3
"""
Created on Wed Apr 17 14:00:26 2024

@author: fabs
TODO: Unite with createPDBs.py
"""

import os
import subprocess
import sys
import urllib.request

import yaml

# where to deploy PDBs
if "NMLDB_DEPLOY_PATH" in os.environ:
    outpath_prefix = os.environ["NMLDB_DEPLOY_PATH"]
else:
    outpath_prefix = os.getcwd()

FAILS = []

if len(sys.argv) > 1:
    files = sys.argv[1:]

    for file in files:
        try:
            path = f"./Databank/Data/Simulations/{file}/"
            outpath = os.path.join(outpath_prefix, file)
            if not os.path.isdir(outpath):
                os.makedirs(outpath, exist_ok=True)
            readme = path + "README.yaml"
            data = yaml.load(open(readme), Loader=yaml.FullLoader)
            if not os.path.isfile(path + data["TPR"][0][0]):
                response = urllib.request.urlretrieve(
                    "https://zenodo.org/record/" + data["DOI"].split(".")[2] + "/files/" + data["TPR"][0][0],
                    path + data["TPR"][0][0],
                )

            if "gromacs" in data["SOFTWARE"]:
                if os.path.isfile(f"{outpath}conf.pdb"):
                    os.remove(f"{outpath}conf.pdb")
                if os.path.isfile(f"{outpath}conf.pdb.gz"):
                    os.remove(f"{outpath}conf.pdb.gz")
                editconf_cmd = ["gmx", "editconf", "-f", f"{path}{data['TPR'][0][0]}", "-o", f"{outpath}conf.pdb"]
                try:
                    subprocess.run(editconf_cmd, check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError as e:
                    raise RuntimeError(f"Command `{' '.join(editconf_cmd)}` failed with stderr: {e.stderr}") from e
                os.remove(path + data["TPR"][0][0])
                gzip_cmd = ["gzip", f"{outpath}conf.pdb"]
                try:
                    subprocess.run(gzip_cmd, check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError as e:
                    raise RuntimeError(f"Command `{' '.join(gzip_cmd)}` failed with stderr: {e.stderr}") from e
        except:
            FAILS.append(file)

    if FAILS:
        print("It was not possible to generate a PDB for:")
        for file in FAILS:
            print(file)

else:
    print("No systems provided")
