"""
`test_analyze` tests functions performing databank recomputing. It starts mostly from
downloading everything from testing repository.

Test data is stored in `./Data/Simulations.1`
"""

from unittest import mock
import os
import glob
import pytest
import logging
import DatabankLib

# Global loaders
# ----------------------------------------------------------------
# we substitute simulation path to test databank analysis scripts
# we load all trajectories from specific zenodo repository
# https://doi.org/10.5281/zenodo.11614468
# On teardown stage we clear the folders.


@pytest.fixture(autouse=True, scope="module")
def header_module_scope():
    _rp = os.path.join(os.path.dirname(__file__), "Data", "Simulations.1")
    with mock.patch.object(DatabankLib, "NMLDB_SIMU_PATH", _rp):
        print("DBG: Mocking simulation path: ", DatabankLib.NMLDB_SIMU_PATH)
        yield
    print("DBG: Mocking completed")


@pytest.fixture(scope="module")
def logger():
    logger = logging.getLogger("test_logger")
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)

    if not logger.handlers:  # Avoid adding multiple handlers during pytest runs
        logger.addHandler(handler)

    yield logger

    # TEARDOWN: clean up handlers after use
    for handler in logger.handlers:
        logger.removeHandler(handler)


@pytest.fixture(scope="module")
def systems():
    from DatabankLib.core import initialize_databank
    s = initialize_databank()
    print(f"Loaded: {len(s)} systems")
    yield s
    # TEARDOWN SYSTEMS
    print('DBG: Wiping temporary calculation data.')
    for _s in s:
        def gbGen(x): return glob.glob(
                    os.path.join(DatabankLib.NMLDB_SIMU_PATH, _s['path'], x))
        clearList = ['*.json', '*.dat', 'conf.gro', 'frame0.gro', '*.dat',
                     '*.buildH', '.*', '#*', '*.def', 'whole.xtc', 'centered.xtc']
        for pat in clearList:
            for f in gbGen(pat):
                os.remove(f)


@pytest.fixture(scope="module")
def systemLoadTraj(systems):
    from DatabankLib.databankLibrary import system2MDanalysisUniverse
    print('DBG: Download trajectory data.')
    for s in systems:
        u = system2MDanalysisUniverse(s)
        del u
    yield
    # TEARDOWN SYSTEM-LOADING
    print('DBG: Wiping trajectory data.')
    for s in systems:
        if s.get('DOI') == 'localhost':
            continue  # don't remove in this case!
        for wp in ['GRO', 'TPR', 'TRJ']:
            try:
                file_ = s[wp][0][0]
            except (KeyError, TypeError):
                continue
            file_ = os.path.join(DatabankLib.NMLDB_SIMU_PATH, s['path'], file_)
            print(file_)
            if os.path.exists(file_):
                os.remove(file_)


# Test functions block.
# ----------------------------------------------------------------
# Every test function is parametrized with system ID to make clear reporting
# about which system actually fails on a test function.


def compareJSONsBtwSD(jsfn: str):
    """Function to call from the test to compare JSON to the precomputed one.
    Raises standard assertion to fail the test.

    Args:
        jsfn (_type_): relative path from simulation folder
    """
    import json
    import numpy as np

    _p1 = os.path.join(os.path.dirname(__file__), "Data", "Simulations.1")
    _p2 = os.path.join(os.path.dirname(__file__), "Data", "Simulations.2")
    jsf1 = os.path.join(_p1, jsfn)
    jsf2 = os.path.join(_p2, jsfn)

    with open(jsf1) as f:
        j1 = json.load(f)
    with open(jsf2) as f:
        j2 = json.load(f)

    assert (type(j1) is type(j2))

    if type(j1) is list:
        assert len(j1) == len(j2), \
            f"Problem in {jsfn} comparison: lists has different lengths!"
        for k1 in range(len(j1)):
            _a = (np.array(j1[k1]) - np.array(j2[k1])).ravel()
            assert (_a**2).sum() < 1e-7, (
               f"Problem in {jsfn} comparison: {k1} line.\n" +
               f"Computed: {str(j1[k1])} \n" +
               f"Pre-computed: {str(j2[k1])}")
    elif type(j1) is dict:
        for k1 in j1:
            _a = (np.array(j1[k1]) - np.array(j2[k1])).ravel()
            assert (_a**2).sum() < 1e-7, (
               f"Problem in {jsfn} comparison: {k1} field.\n" +
               f"Computed: {str(j1[k1])} \n" +
               f"Pre-computed: {str(j2[k1])}")

    print(f"Data {jsfn} was compared against precomputed!")


@pytest.mark.parametrize("systemid", [86, 243, 281, 566, 787])
def test_analyze_apl(systems, systemLoadTraj, systemid, logger):
    from DatabankLib.analyze import computeAPL

    s = systems.loc(systemid)

    rCode = computeAPL(s, logger)
    assert rCode == DatabankLib.analyze.RCODE_COMPUTED

    cFile = os.path.join(DatabankLib.NMLDB_SIMU_PATH,
                         s['path'], 'apl.json')
    assert os.path.isfile(cFile)
    assert os.path.getsize(cFile) > 1e3
    compareJSONsBtwSD(
        os.path.relpath(cFile, DatabankLib.NMLDB_SIMU_PATH)
    )


@pytest.mark.parametrize("systemid, rcodex",
                         [(281, DatabankLib.RCODE_COMPUTED),
                          (566, DatabankLib.RCODE_COMPUTED),
                          (787, DatabankLib.RCODE_ERROR),
                          (243, DatabankLib.RCODE_COMPUTED),
                          (86,  DatabankLib.RCODE_COMPUTED)])
def test_analyze_op(systems, systemLoadTraj, systemid, rcodex, logger):
    from DatabankLib.analyze import computeOP
    from DatabankLib.settings.molecules import lipids_dict
    s = systems.loc(systemid)
    rCode = computeOP(s, logger)
    assert rCode == rcodex
    if rcodex == DatabankLib.RCODE_ERROR:
        return

    for lip in set(s['COMPOSITION'].keys()).intersection(set(lipids_dict.keys())):
        cFile = os.path.join(DatabankLib.NMLDB_SIMU_PATH,
                             s['path'], lip + 'OrderParameters.json')
        assert os.path.isfile(cFile), f"File {cFile} wasn't created for {lip}!"
        assert os.path.getsize(cFile) > 1e3, f"File {cFile} for {lip} is less than 1K!"
        compareJSONsBtwSD(
            os.path.relpath(cFile, DatabankLib.NMLDB_SIMU_PATH)
        )


@pytest.mark.parametrize("systemid, rcodex",
                         [(281, DatabankLib.RCODE_COMPUTED),
                          (566, DatabankLib.RCODE_COMPUTED),
                          (787, DatabankLib.RCODE_ERROR),
                          (243, DatabankLib.RCODE_COMPUTED),
                          (86,  DatabankLib.RCODE_COMPUTED)],)
def test_analyze_ff(systems, systemLoadTraj, systemid, rcodex, logger):
    from DatabankLib.analyze import computeFF
    s = systems.loc(systemid)
    rCode = computeFF(s, logger)
    assert rCode == rcodex
    if rcodex == DatabankLib.RCODE_ERROR:
        return

    for fn in ['FormFactor.json', 'TotalDensity.json', 'WaterDensity.json',
               'LipidDensity.json']:
        cFile = os.path.join(DatabankLib.NMLDB_SIMU_PATH,
                             s['path'], fn)
        assert os.path.isfile(cFile)
        assert os.path.getsize(cFile) > 1e3
        compareJSONsBtwSD(
            os.path.relpath(cFile, DatabankLib.NMLDB_SIMU_PATH)
        )


@pytest.mark.parametrize("systemid, rcodex",
                         [(281, DatabankLib.RCODE_COMPUTED),
                          (566, DatabankLib.RCODE_COMPUTED),
                          (787, DatabankLib.RCODE_ERROR),
                          (243, DatabankLib.RCODE_COMPUTED),
                          (86,  DatabankLib.RCODE_COMPUTED)])
def test_analyze_nmrpca(systems, systemLoadTraj, systemid, rcodex, logger):
    from DatabankLib.analyze import computeNMRPCA
    s = systems.loc(systemid)
    rCode = computeNMRPCA(s, logger)
    assert rCode == rcodex

    if rCode == DatabankLib.RCODE_ERROR:
        return

    cFile = os.path.join(DatabankLib.NMLDB_SIMU_PATH,
                         s['path'], 'eq_times.json')
    assert os.path.isfile(cFile)
    assert os.path.getsize(cFile) > 10
    compareJSONsBtwSD(
        os.path.relpath(cFile, DatabankLib.NMLDB_SIMU_PATH)
    )
