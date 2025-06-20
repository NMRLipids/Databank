# Test the Molecules class
import pytest
import os



# run only on sim2 mocking data
pytestmark = pytest.mark.sim2

def test_initialize_n():
    import DatabankLib
    from DatabankLib.core import initialize_databank
    from DatabankLib.settings.molecules import lipids_set
    #from molecules import lipids, nonlipids,lipids_set, nonlipids_set
    print(DatabankLib.NMLDB_MOL_PATH)
    print(lipids_set)

    assert len(lipids_set) == 5, "LipidSet should have length 5"
    # Access element 'POPE' in lipids_set
    assert 'POPE' in lipids_set, "POPE should be in lipids_set"
    pope = lipids_set.get('POPE')
    assert pope is not None, "POPE should not be None"
    assert pope.name == 'POPE', "POPE name should be 'POPE'" 
    print(f"POPE metadata: {pope._metadata}")   
    
