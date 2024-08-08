from unittest import mock
import os, glob
import pytest
import DatabankLib

## Global loaders
## ----------------------------------------------------------------
## we substitute simulation path to test databank analysis scripts
## we load all trajectories from specific zenodo repository
## https://doi.org/10.5281/zenodo.11614468
## On teardown stage we clear the folders.


@pytest.fixture(autouse=True, scope="session")
def header_session_scope():
    _rp = os.path.join(os.path.dirname(__file__), "Data", "Simulations.1")
    with mock.patch.object(DatabankLib, "NMLDB_SIMU_PATH", _rp):
        print("DBG: Mocking simulation path: ", DatabankLib.NMLDB_SIMU_PATH)
        yield
    print("DBG: Mocking completed")

@pytest.fixture(scope="module")
def systems():
    from DatabankLib.core import initialize_databank
    s = initialize_databank()
    print(f"Loaded: {len(s)} systems")
    yield s
    ## TEARDOWN SYSTEMS
    print('DBG: Wiping temporary calculation data.')
    for _s in s:
        gbGen = lambda x: glob.glob(os.path.join(DatabankLib.NMLDB_SIMU_PATH, _s['path'], x))
        clearList = ['*.json', '*.dat', 'conf.gro', 'frame0.gro', '*.dat', '*.buildH', '.*', '#*',
                     '*.def', 'whole.xtc']
        for pat in clearList:
            for f in gbGen(pat):
                os.remove(f)


@pytest.fixture(scope="module")
def systemLoadTraj(systems):
    from DatabankLib.databankLibrary import system2MDanalysisUniverse
    print('DBG: Download trajectory data.')
    for s in systems:
        u = system2MDanalysisUniverse(s)
    yield
    ## TEARDOWN SYSTEM-LOADING
    print('DBG: Wiping trajectory data.')
    for s in systems:
        for wp in ['GRO', 'TPR', 'TRJ']:
            try:
                file_ = s[wp][0][0]
            except (KeyError, TypeError):
                continue
            file_ = os.path.join(DatabankLib.NMLDB_SIMU_PATH, s['path'], file_)
            print(file_)
            if os.path.exists(file_):
                os.remove(file_)

## Test functions block.
## ----------------------------------------------------------------
## Every test function is parametrized with system ID to make clear reporting
## about which system actually fails on a test function.

@pytest.mark.parametrize("systemid", [281, 566, 787])
def test_analyze_apl(systems, systemLoadTraj, systemid):
    from DatabankLib.analyze import computeAPL

    s = systems.loc(systemid)

    rCode = computeAPL(s)
    assert rCode == DatabankLib.analyze.RCODE_COMPUTED

    cFile = os.path.join(DatabankLib.NMLDB_SIMU_PATH, 
                                    s['path'], 'apl.json')
    assert os.path.isfile(cFile)
    assert os.path.getsize(cFile) > 1e3


@pytest.mark.parametrize("systemid, rcodex", 
                         [   (281, DatabankLib.RCODE_COMPUTED), 
                             (566, DatabankLib.RCODE_COMPUTED), 
                             (787, DatabankLib.RCODE_ERROR),
                             (243, DatabankLib.RCODE_COMPUTED),
                             (86,  DatabankLib.RCODE_COMPUTED)     ], )
def test_analyze_op(systems, systemLoadTraj, systemid, rcodex):
    from DatabankLib.analyze import computeOP
    s = systems.loc(systemid)
    rCode = computeOP(s)
    assert rCode == rcodex
    if rcodex == DatabankLib.RCODE_ERROR:
        return
    
    # now check only 1st lipid
    cFile = os.path.join(DatabankLib.NMLDB_SIMU_PATH, 
        s['path'], list(s['COMPOSITION'].keys())[0] + 'OrderParameters.json')
    assert os.path.isfile(cFile)
    assert os.path.getsize(cFile) > 1e3 

