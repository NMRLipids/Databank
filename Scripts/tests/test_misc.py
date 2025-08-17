"""
`test_misc` contains unit tests of auxiliary functions.
Test data is stored in `./Data/Simulations.2`

-------------------------------------------------------------------------------
NOTE: globally import of DatabankLib is **STRICTLY FORBIDDEN** because it 
      breaks the substitution of global path folders
"""

import os
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


def test_maicos_interface(tmpdir):
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader
    from DatabankLib.maicos import DensityPlanar
    import numpy as np

    folder = str(tmpdir)
    # Skip unwrap/pack for speed - trajectories are already centered and whole
    base_options = {"unwrap": False, "bin_width": 1, "pack": False}
    zlim = {"zmin": 0, "zmax": 8.67623}
    dens_options = {**zlim, **base_options}

    pdb_content = \
"""TITLE     Single 4-site water
CRYST1   63.646   63.646   86.762  90.00  90.00  90.00 P 1           1
ATOM      1  OW  wate    1       4.370  31.010  61.030  1.00  0.00            
ATOM      2  HW1 wate    1       3.800  30.820  61.780  1.00  0.00            
ATOM      3  HW2 wate    1       5.250  30.760  61.320  1.00  0.00            
ATOM      4  MW  wate    1       4.410  30.950  61.170  1.00  0.00  
END
"""
    tmpfile = os.path.join(folder, "test_water_maicos.pdb")
    with open(tmpfile, 'w') as f:
        f.write(pdb_content)
    u = mda.Universe(
        tmpfile,
        format="PDB",
    )
    nul_dim = u.dimensions
    # Store coordinates for each frame
    n_frames = 10
    trajectory = np.zeros((n_frames, u.atoms.n_atoms, 3), dtype=np.float32)

    for i in range(n_frames):
        coords = u.atoms.positions.copy()
        # move all atoms a little bit in x direction
        coords = np.add(coords, [0.1 * i, 0.0, 0.0])  
        trajectory[i] = coords

    # Assign trajectory to Universe
    u.load_new(trajectory, order='fac', format=MemoryReader)
    for ts in u.trajectory:
        ts.dimensions = nul_dim
    # Now we are done!
    u.add_TopologyAttr("elements")
    u.atoms.elements = ['O', 'H', 'H', 'X']  # Assign elements manually

    # Simulate MAICoS calls
    cos_water = DensityPlanar(
        u.atoms,
        dens="electron",
        **dens_options,
        output=os.path.join(folder, "DiporderWater.json"),
    )
    cos_water.run()
    assert True
