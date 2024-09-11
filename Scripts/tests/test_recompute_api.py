"""
`test_recompute_api.py` tests those parts of API, which require building MDAnalysis universe 
and from-trajectories computations. 

Test folder: Data/Simulations.1
"""

from unittest import mock
import os, glob
from tempfile import TemporaryDirectory
from contextlib import nullcontext as does_not_raise
import pytest

## Global fixtures
## ----------------------------------------------------------------
## we substitute simulation path to test databank analysis scripts

@pytest.fixture(autouse=True, scope="module")
def header_module_scope():
    import DatabankLib
    _rp = os.path.join(os.path.dirname(__file__), "Data", "Simulations.1")
    with mock.patch.object(DatabankLib, "NMLDB_SIMU_PATH", _rp) as dbl:
        print("DBG: Mocking simulation path: ", DatabankLib.NMLDB_SIMU_PATH)
        yield
    print("DBG: Mocking completed")

@pytest.fixture(scope="module")
def systems():
    import DatabankLib
    from DatabankLib.core import initialize_databank
    s = initialize_databank()
    print(f"Loaded: {len(s)} systems")
    yield s
    ## TEARDOWN SYSTEMS
    print('DBG: Wiping temporary calculation data.')
    for _sid in [787]:
        _s = s.loc(787)
        gbGen = lambda x: glob.glob(os.path.join(DatabankLib.NMLDB_SIMU_PATH, _s['path'], x))
        clearList = ['*.xtc', '*.gro']
        for pat in clearList:
            for f in gbGen(pat):
                os.remove(f)


## Test functions block.
## ----------------------------------------------------------------
## Every test function is parametrized with system ID to make clear reporting
## about which system actually fails in a test function.

@pytest.mark.parametrize("systemid, natoms, nframes", 
                         [ (243, 12208, 121),   # with TPR, united-atom, local
                           (787, 46290, 1141),  # with GRO only, aa, network
                         ]  )
def test_system2MDAnalysisUniverse(systems, systemid, natoms, nframes):
    from DatabankLib.databankLibrary import system2MDanalysisUniverse
    s = systems.loc(systemid)
    u = system2MDanalysisUniverse(s)
    assert len(u.atoms) == natoms
    assert u.trajectory.n_frames == nframes
    with does_not_raise() as e:
        # check that it doesn't fail iterating over the trajectory
        for f in u.trajectory:
            pass

@pytest.fixture(scope='function')
def failSys():
    import DatabankLib
    with TemporaryDirectory(prefix=DatabankLib.NMLDB_SIMU_PATH + os.sep) as tmpd:
        s = {
            'DOI': 'localhost',
            'GRO': [['md.gro']],
            'TRJ': [['md.trr']],
            'path': os.path.relpath(DatabankLib.NMLDB_SIMU_PATH, tmpd),
            'SOFTWARE': 'gromacs'
        }
        yield s
    # TEARDOWN
    pass

@pytest.mark.xfail(reason="Localhost with non-downloaded files",
                   raises=FileNotFoundError)
def test_fail1_system2MDAnalysisUniverse(failSys):
    from DatabankLib.databankLibrary import system2MDanalysisUniverse
    u = system2MDanalysisUniverse(failSys)

def hashFV(x):
    import numpy as np
    a = np.array(x)
    a = np.around(a, 4)
    a *= 1e4
    a = np.array(a, dtype='int32')
    a = tuple(a.tolist())
    return hash(a)

@pytest.mark.parametrize("systemid, lipid, fvhash", 
                         [ (243, "DPPC", -5227956720741036084),  # with TPR, united-atom, local
                           (787, "POPC", 4799549858726566566)    # with GRO only, aa, network
                         ]  )
def test_PJangle(systems, systemid, lipid, fvhash):
    from DatabankLib.databankLibrary import (read_trj_PN_angles, 
        system2MDanalysisUniverse, simulation2universal_atomnames)
    s = systems.loc(systemid)
    u = system2MDanalysisUniverse(s)
    a1 = simulation2universal_atomnames(s, lipid, "M_G3P2_M")
    a2 = simulation2universal_atomnames(s, lipid, "M_G3N6_M")

    a,b,c,d = read_trj_PN_angles(lipid, a1, a2, u)
    # time-molecule arrays
    assert len(a) == sum(s['COMPOSITION'][lipid]['COUNT'])
    # time-averaged list
    assert len(b) == sum(s['COMPOSITION'][lipid]['COUNT'])
    # comparing per-molecule hash. Suppose that's enough
    assert hashFV(b) == fvhash
    # overall mean
    assert c > 0
    # overall std
    assert d > 0 
