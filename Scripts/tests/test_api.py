from unittest import mock
import sys, os
from urllib.error import HTTPError
import pytest

import DatabankLib

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

def test_initialize_n(systems):
    assert len(systems) == 4

def test_print_README(systems, capsys):
    from DatabankLib.core import print_README
    sys0 = systems[0]
    print_README(sys0)
    output: str = capsys.readouterr().out.rstrip()
    fDOI = output.find('DOI:') != -1
    fTEMP = output.find('TEMPERATURE:') != -1
    assert fDOI and fTEMP

def test_CalcAreaPerMolecule(systems):
    from DatabankLib.databankLibrary import CalcAreaPerMolecule
    sys0 = systems[0]
    apm = CalcAreaPerMolecule(sys0)
    print(apm) 
    assert 1
