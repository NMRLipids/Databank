#!/usr/bin/env python3
# coding: utf-8

"""
:program: makeRanking.py
:description: creates different types of ranking lists of simulations based on their
              quality against experiments. The ranking lists are stored in
              ``Data/Ranking/`` folder in JSON format. The lists can be shown
              with the ``plotQuality.ipynb``
"""

import os
import json

from DatabankLib import NMLDB_DATA_PATH, NMLDB_SIMU_PATH
from DatabankLib.databankLibrary import lipids_dict
from DatabankLib.core import initialize_databank
from DatabankLib.jsonEncoders import CompactJSONEncoder

if __name__ == "__main__":
    systems = initialize_databank()

    # ---- Making list of qualities
    qualities = []
    for system in systems:
        quality_dict = {}
        path = os.path.join(NMLDB_SIMU_PATH, system['path'])

        total_quality_file_path = os.path.join(path, 'SYSTEM_quality.json')

        fragment_Q = {}
        for lipid in system['COMPOSITION']:
            quality_file = os.path.join(path, lipid + '_FragmentQuality.json')
            try:
                with open(quality_file, 'r') as json_file:
                    fragment_Q[lipid] = json.load(json_file)
            except Exception:
                continue

        if (os.path.isfile(total_quality_file_path)):
            with open(total_quality_file_path, 'r') as json_file:
                fragment_Q['TotalQuality'] = json.load(json_file)

        ff_quality_fpath = os.path.join(path, 'FormFactorQuality.json')
        if (os.path.isfile(ff_quality_fpath)):
            with open(ff_quality_fpath) as json_file:
                FFq = json.load(json_file)
            try:
                fragment_Q['TotalQuality']
            except KeyError:
                fragment_Q['TotalQuality'] = {}
            try:
                fragment_Q['TotalQuality']['FFQuality'] = FFq[0]
            except (KeyError, TypeError):
                fragment_Q['TotalQuality']['FFQuality'] = FFq
            json_file.close()

        fragment_Q['system'] = system
        qualities.append(fragment_Q)

    # ---- Sort based on total quality of a simulation
    fragments = ['total', 'tails', 'headgroup']

    for sort_based_on in fragments:
        new_qualities = []
        for i in qualities:
            try:
                if i['TotalQuality'][sort_based_on] > 0:
                    new_qualities.append(i)
            except (KeyError, TypeError):
                continue

        sorted_qualities = sorted(
            new_qualities, key=lambda i: i['TotalQuality'][sort_based_on], reverse=True)

        outputfile = os.path.join(NMLDB_DATA_PATH, 'Ranking',
                                  'SYSTEM_' + sort_based_on + '_Ranking.json')
        with open(outputfile, "w") as fp:
            json.dump(sorted_qualities, fp, default=str, cls=CompactJSONEncoder)
        print(f"Sorted based on {sort_based_on} quality and saved to {outputfile}")

    new_qualities = []
    for i in qualities:
        try:
            if i['TotalQuality']['FFQuality'] > 0:
                new_qualities.append(i)
        except (KeyError, TypeError):
            continue

    sorted_qualities = sorted(new_qualities,
                              key=lambda i: i['TotalQuality']['FFQuality'])

    outputfile = os.path.join(NMLDB_DATA_PATH, 'Ranking',
                              'SYSTEM_FormFactor_Ranking.json')
    with open(outputfile, "w") as fp:
        json.dump(sorted_qualities, fp, default=str, cls=CompactJSONEncoder)
    print('Sorted based on form factor quality and saved to', outputfile)

    # ---- Sorting best simulations for each lipid
    fragments = ['total', 'sn-1', 'sn-2', 'headgroup']

    for sort_based_on in fragments:
        for lipid in lipids_dict:
            new_qualities = []
            for i in qualities:
                try:
                    if i[lipid][sort_based_on] > 0:
                        new_qualities.append(i)
                except (KeyError, TypeError):
                    continue

            sorted_qualities = sorted(
                new_qualities, key=lambda i: i[lipid][sort_based_on], reverse=True)

            if sorted_qualities:
                outputfile = os.path.join(NMLDB_DATA_PATH, 'Ranking',
                                          lipid + '_' + sort_based_on + '_Ranking.json')
                with open(outputfile, "w") as fp:
                    json.dump(sorted_qualities, fp, default=str, cls=CompactJSONEncoder)
                print(f"Quality of {sort_based_on} of {lipid} "
                      f"sorted and saved to {outputfile}")
