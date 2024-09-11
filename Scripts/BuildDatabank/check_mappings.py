#!/usr/bin/env python3
"""
It checks mapping files, but going only through systems
which have TPR in folders preloaded. That is required for
debugging mappings code.
"""

import os
import sys
import MDAnalysis as mda
from tqdm import tqdm
import traceback

sys.path.append('..')
import DatabankLib
from DatabankLib import NMLDB_SIMU_PATH
from DatabankLib.core import initialize_databank
from DatabankLib.databankLibrary import lipids_dict, loadMappingFile

if __name__ == "__main__":
    systems = initialize_databank()

    for system in tqdm(systems):
        try:
            if 'TPR' not in system.keys():
                sys.stderr.write(f"Skipping {system['SYSTEM']} because there is no TPR\n")
                continue
            tpr_name = os.path.join(NMLDB_SIMU_PATH, system['path'], system['TPR'][0][0])
            if not os.path.exists(tpr_name):
                sys.stderr.write(f"""
        Skipping {system["SYSTEM"]} because TPR is not downloaded
        ==> {tpr_name}
        """)
                continue
            print(f"Checking system {system['path']} ...")
            # opening Universe using TPR only
            try:
                tpars = mda.topology.TPRParser.TPRParser(tpr_name)
                tpl = tpars.parse()
                u = mda.Universe(tpl)
            except Exception as e:
                sys.stderr.write(f"Problem openening {tpr_name} with MDAnalysis!\n")
                sys.stderr.write(str(e))
                continue
            for molecule in system['COMPOSITION']:
                errnum = 0
                m_file = system['COMPOSITION'][molecule]['MAPPING']
                mapping_dict = loadMappingFile(m_file)
                # go over all records
                for mk in mapping_dict.keys():
                        selection = 'resname ' + system['COMPOSITION'][molecule]['NAME'] + ' and name ' +  mapping_dict[mk]['ATOMNAME']
                        NatomsFromMapping = len(u.select_atoms(selection))
                        NatomsFromREADME = system['COMPOSITION'][molecule]['COUNT']
                        if NatomsFromMapping !=  NatomsFromREADME and molecule not in lipids_dict:
                            print(f"""
        Found problematic system: {system['SYSTEM']}
        Molecule named {system['COMPOSITION'][molecule]['NAME']}
        MDA selection: "{selection}"
        Atoms from mapping/in readme: {NatomsFromMapping} / {NatomsFromREADME}
        """)
                            errnum += 1
                # end mapping-dict cycle
                if errnum == 0:
                    print("Everything is OK for molecule %s!" % molecule)
        except Exception as e:
            sys.stderr.write(f"Unexpected error in system {system['ID']}\n")
            sys.stderr.write(traceback.format_exc())
