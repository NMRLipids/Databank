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
    print(f"POPE metadata: {pope.metadata}")  
    metadata = pope.metadata
    assert metadata is not None, "POPE metadata should not be None"
    assert metadata['bioschema_properties'] is not None, "POPE bioschema_properties should not be None"
    assert metadata['NMRlipids']['id'] == 'POPE', "POPE bioschema_properties name should be 'POPE'"
    assert metadata['NMRlipids']['name'] == '1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine', "POPE NMRlipids name should match expected value"
    assert metadata['bioschema_properties']['molecularWeight'] == 718, "POPE bioschema_properties molecularWeight should be 718.0"  
    assert metadata['bioschema_properties']['molecularFormula'] == 'C39H76NO8P', "POPE bioschema_properties formula should be 'C39H76NO8P'"
    assert metadata['bioschema_properties']['smiles'] == 'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OCCN)OC(=O)CCCCCCC/C=C\CCCCCCCC', "POPE bioschema_properties smiles should match expected value"