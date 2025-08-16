"""
`test_misc` contains unit tests of auxiliary functions.
Test data is stored in `./Data/Simulations.2`

-------------------------------------------------------------------------------
NOTE: globally import of DatabankLib is **STRICTLY FORBIDDEN** because it 
      breaks the substitution of global path folders
"""

import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = pytest.mark.nodata


def test_uname2element():
    """Test uname2element function."""
    from DatabankLib.settings.elements import uname2element

    check.equal(uname2element("M_C1_M"), "C")
    check.equal(uname2element("M_G1_M"), "C")
    check.equal(uname2element("M_C1N3_M"), "N")
    check.equal(uname2element("M_X_M"), "X")

    with pytest.raises(KeyError):
        uname2element("UnknownElement")
