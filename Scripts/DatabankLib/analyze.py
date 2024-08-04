
import json
import os
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)

from . import *
from .databank_defs import *
from .core import *
from .databankLibrary import GetNlipids, system2MDanalysisUniverse
from .jsonEncoders import CompactJSONEncoder

#TODO: use in calcAPL.py
def computeAPL(system: dict, recompute: bool = False):
    """Generate apl.json analysis file for a system.

    Args:
        system (dict): one of systems of the Databank
        recompute (bool, optional): Delete previous apl.json and recompute it if True. Defaults to False.
    """
    ## reading software and file path for simulations
    software = system['SOFTWARE']
    path = system['path']

    ## this is checking if area per lipid is already calculated for the systems
    outfilename = os.path.join(NMLDB_SIMU_PATH,  path, 'apl.json')
    if os.path.isfile(outfilename):
        if recompute:
            os.unlink(outfilename)
        else:
            return
    
    print('Analyzing: ', path)
    print('Will write into: ', outfilename)

    ## calculates the total number of lipids
    Nlipid = GetNlipids(system)

    ## makes MDAnalysis universe from the system. This also downloads the data if not yet locally available
    u = system2MDanalysisUniverse(system)

    if u is None:
        print('Generation of MDAnalysis universe failed in folder', path)
        return
    
    ## this calculates the area per lipid as a function of time and stores it in the databank
    apl = {}
    for ts in tqdm(u.trajectory, desc='Scanning the trajectory'):
        if u.trajectory.time >= system['TIMELEFTOUT']*1000:
            dimensions = u.dimensions
            aplFrame = u.dimensions[0]*u.dimensions[1]*2/Nlipid
            apl[u.trajectory.time] = aplFrame

    with open(outfilename, 'w') as f:
        json.dump(apl, f, cls=CompactJSONEncoder)
    
    return

