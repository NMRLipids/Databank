from unittest import mock
import sys, os
from urllib.error import HTTPError
import pytest

import DatabankLib

@pytest.fixture(autouse=True, scope="session")
def header_session_scope():
    _rp = os.path.join(os.path.dirname(__file__), "Data", "Simulations.2")
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

def test_initialize_n(systems):
    assert len(systems) == 5

def test_print_README(systems, capsys):
    from DatabankLib.core import print_README
    sys0 = systems[0]
    print_README(sys0)
    output: str = capsys.readouterr().out.rstrip()
    fDOI = output.find('DOI:') != -1
    fTEMP = output.find('TEMPERATURE:') != -1
    assert fDOI and fTEMP

@pytest.mark.parametrize("systemid, result", 
                         [   (281, 64.722), 
                             (566, 61.306), 
                             (787, 78.237),
                             (243, 62.276),
                             (86,  60.460)     ], )
def test_CalcAreaPerMolecule(systems, systemid, result):
    from DatabankLib.databankLibrary import CalcAreaPerMolecule
    sys0 = systems.loc(systemid)
    apm = CalcAreaPerMolecule(sys0)
    assert abs(apm - result) < 6e-4


@pytest.mark.parametrize("systemid, result", 
                         [   (281, 128), 
                             (566, 128), 
                             (787, 120),
                             (243, 72),
                             (86,  128)     ], )
def test_GetNLipids(systems, systemid, result):
    from DatabankLib.databankLibrary import GetNlipids
    sys0 = systems.loc(systemid)
    nlip = GetNlipids(sys0)
    assert nlip == result


@pytest.mark.parametrize("systemid, result", 
                         [   (281, [0.264, 0.415, 0.614]), 
                             (566, [0.263, 0.404, 0.604]), 
                             (243, [0.281, 0.423, 0.638]),
                             (86,  [0.266, 0.421, 0.63])     ], )
def test_GetFormFactorMin(systems, systemid, result):
    from DatabankLib.databankLibrary import GetFormFactorMin
    import numpy as np
    sys0 = systems.loc(systemid)
    ffl = GetFormFactorMin(sys0)
    err = ( (np.array(ffl[:3]) - np.array(result))**2 ).sum()
    assert err < 1e-9

@pytest.mark.parametrize("systemid, result", 
                         [   (281, 31.5625), 
                             (566, 31.0), 
                             (787, 75.0),
                             (243, 39.7778),
                             (86,  27.75)     ], )
def test_getHydrationLevel(systems, systemid, result):
    from DatabankLib.databankLibrary import getHydrationLevel
    sys0 = systems.loc(systemid)
    hl = getHydrationLevel(sys0)
    assert abs(hl - result) < 1e-4


@pytest.mark.parametrize("systemid, lipid, result", 
                         [   (281, ['POPC'], [1]), 
                             (566, ['POPC','CHOL'], [.9375,.0625]), 
                             (787, ['TOCL', 'POPC', 'POPE'], [0.25,0.5,0.25]),
                             (243, ['DPPC'], [1]),
                             (86,  ['POPE'], [1])     ], )
def test_calcLipidFraction(systems, systemid, lipid, result):
    from DatabankLib.databankLibrary import calcLipidFraction
    sys0 = systems.loc(systemid)
    assert calcLipidFraction(sys0, 'SOPC') == 0 # absent lipid
    err = 0
    i = 0
    for lip in lipid:
        err += ( calcLipidFraction(sys0, lip) - result[i] )**2
        i += 1
    assert err < 1e-4

@pytest.mark.parametrize("systemid, result", 
                         [   (281, [-0.1610, -0.1217] ), 
                             (566, [-0.1714, -0.1142]), 
                             (243, [-0.1764, -0.1784]),
                             (86,  [-0.1933, -0.1568])     ], )
def test_averageOrderParameters(systems, systemid, result):
    from DatabankLib.databankLibrary import averageOrderParameters
    sys0 = systems.loc(systemid)
    sn1,sn2 =  averageOrderParameters(sys0)
    assert (sn1-result[0])**2 + (sn2-result[1])**2 < 1e-5


@pytest.mark.parametrize("systemid", [ 787], )
def test_raises_averageOrderParameters(systems, systemid):
    from DatabankLib.databankLibrary import averageOrderParameters
    sys0 = systems.loc(systemid)
    with pytest.raises(FileNotFoundError) as exc_info:
        sn1,sn2 =  averageOrderParameters(sys0)
    assert 'OrderParameters.json' in str(exc_info.value)
