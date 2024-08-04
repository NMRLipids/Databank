from unittest import mock
import os
import pytest
import DatabankLib

@pytest.fixture(autouse=True, scope="session")
def header_session_scope():
    _rp = os.path.join(os.path.dirname(__file__), "Data", "Simulations")
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

def test_analyze_apl(systems):
    from DatabankLib.analyze import computeAPL
    aplCreated = []
    aplNZ = []
    for s in systems:
        computeAPL(s)
        cFile = os.path.join(DatabankLib.NMLDB_SIMU_PATH, 
                                        s['path'], 'apl.json')
        aplCreated.append( os.path.isfile(cFile) )
        aplNZ.append( os.path.getsize(cFile) > 1e3 )
    assert all(aplCreated) and all(aplNZ)