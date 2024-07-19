"""
Databank routines designed for calling from inside jupyter notebooks. 
"""

import pandas as pd
import numpy as np
from IPython.display import display
from IPython.display import Markdown

def ShowTable(SortedQualities, quality):
    """
    Shows a table of simulation qualities against experimental data.

    :param SortedQualities: list of dictionaries to be shown, available in folder ``Data/Ranking/``
    :param quality: should be either ``TotalQuality`` or universal lipid name. First one shows the system total quality. Latter shows the individual lipid quality.

    """
    rounding = ['headgroup', 'sn-1', 'sn-2', 'total', 'tails', 'FFQuality']
    QualityTable = []
    pd.set_option('display.max_rows', None)
    for i in SortedQualities:
        StoredToTable = []
        for k, v in i[quality].items():
            if k in rounding:
                if v and v != float("inf") and not np.isnan(v):
                    i[quality][k] = round(float(v), 2)

        StoredToTable = i[quality]
        try:
            StoredToTable['Forcefield'] = i['system']['FF']
        except:
            display(Markdown('**FAILURE:** no FF defined for the system'))
            display(i['system']['path'])
            continue

        molecules = ''
        MolNumbers = ''
        for lipid in i['system']['COMPOSITION']:
            molecules = molecules + lipid + ':'
            MolNumbers = MolNumbers + str(np.sum(i['system']['COMPOSITION'][lipid]['COUNT']))  + ':'
        StoredToTable['Molecules'] = molecules[:-1]
        StoredToTable['Number of molecules'] = ' (' + MolNumbers[:-1] + ')'
        StoredToTable['Temperature'] = i['system']['TEMPERATURE']
        StoredToTable['ID'] = i['system']['ID']
        QualityTable.append(StoredToTable)    
    display(pd.json_normalize(QualityTable))

