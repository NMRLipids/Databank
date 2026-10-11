"""Unit tests for MAICoS-specific functionality."""

import json
import logging
import os
import shutil
import tempfile

import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = [pytest.mark.sim1, pytest.mark.min]


@pytest.fixture
def tip4p_trajectory(tmpdir):
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.coordinates.memory import MemoryReader

    folder = str(tmpdir)
    pdb_content = """TITLE     Single 4-site water
CRYST1   63.646   63.646   86.762  90.00  90.00  90.00 P 1           1
ATOM      1  OW  wate    1       4.370  31.010  61.030  1.00  0.00
ATOM      2  HW1 wate    1       3.800  30.820  61.780  1.00  0.00
ATOM      3  HW2 wate    1       5.250  30.760  61.320  1.00  0.00
ATOM      4  MW  wate    1       4.410  30.950  61.170  1.00  0.00
END
"""
    tmpfile = os.path.join(folder, "test_water_maicos.pdb")
    with open(tmpfile, "w") as f:
        f.write(pdb_content)
    with pytest.warns(UserWarning, match="Element information is missing"):
        u = mda.Universe(tmpfile, format="PDB")
    nul_dim = u.dimensions
    trajectory = np.zeros((10, u.atoms.n_atoms, 3), dtype=np.float32)
    shift_v = np.array([0.0, 0.0, -0.5], dtype=np.float32)
    for i in range(10):
        trajectory[i] = u.atoms.positions.copy() + shift_v * i
    u.load_new(trajectory, order="fac", format=MemoryReader)
    for ts in u.trajectory:
        ts.dimensions = nul_dim
    yield u, folder
    shutil.rmtree(folder)


def test_maicos_interface(tip4p_trajectory, tmp_path):
    from fairmd.lipids.analib.maicos import DensityPlanar

    u, _ = tip4p_trajectory
    u.add_TopologyAttr("elements")
    u.atoms.elements = ["O", "H", "H", "Dummy"]
    ofname = tmp_path / "output.json"
    dp = DensityPlanar(
        u.atoms,
        dens="electron",
        zmin=0,
        zmax=20,
        unwrap=False,
        bin_width=1,
        pack=False,
        output=ofname,
    )
    dp.run()
    dp.save()
    assert os.path.isfile(ofname)

    with open(ofname) as output_file:
        output = json.load(output_file)
    assert len(output) == 20  # noqa: PLR2004
    assert all(len(row) == 3 for row in output)  # noqa: PLR2004
    assert max(row[1] for row in output) > 0
    check.equal(output[0][1], pytest.approx(0.0, abs=1e-6))
    check.equal(output[14][1], pytest.approx(0.493686, abs=1e-6))
    check.equal(output[14][2], pytest.approx(0.294148, abs=1e-6))
    check.equal(output[-1][1], pytest.approx(0.0, abs=1e-6))


@pytest.mark.parametrize(
    ("system", "expected"),
    [
        ({"TYPEOFSYSTEM": "miscellaneous"}, False),
        ({"TYPEOFSYSTEM": "lipid bilayer"}, True),
        ({"TYPEOFSYSTEM": "lipid bilayer", "WARNINGS": {"ORIENTATION": "warning"}}, False),
        ({"TYPEOFSYSTEM": "lipid bilayer", "WARNINGS": {"PBC": "hexagonal-box"}}, False),
        ({"TYPEOFSYSTEM": "lipid bilayer", "WARNINGS": {"NOWATER": True}}, False),
        ({"TYPEOFSYSTEM": "lipid bilayer", "WARNINGS": {"PBC": "rectangular-box"}}, True),
        ({"TYPEOFSYSTEM": "lipid bilayer", "WARNINGS": {}}, True),
    ],
)
def test_maicos_suitable(system, expected):
    from fairmd.lipids.analib.maicos import is_system_suitable_4_maicos

    assert is_system_suitable_4_maicos(system) is expected


def test_maicos_what_to_compute(caplog, logger):
    from fairmd.lipids import FMDL_SIMU_PATH, RCODE_ERROR
    from fairmd.lipids.analyze import computeMAICOS
    from fairmd.lipids.core import System

    with tempfile.TemporaryDirectory(dir=FMDL_SIMU_PATH) as path:
        for filename in ("file.xtc", "file.tpr", "file.pdb"):
            open(os.path.join(path, filename), "w").close()

        system = System(
            data={
                "DOI": "00.0000/abcd",
                "path": os.path.basename(path),
                "TYPEOFSYSTEM": "lipid bilayer",
                "SOFTWARE": "GROMACS",
                "TPR": [["file.tpr"]],
                "TRJ": [["file.xtc"]],
                "ID": 9999,
                "COMPOSITION": {
                    "SOL": {"NAME": "SPC", "MAPPING": "mappingSPCwater.yaml", "COUNT": 1000},
                    "DPPC": {"NAME": "DPPC", "MAPPING": "mappingDPPCberger.yaml", "COUNT": 100},
                },
            },
        )
        with caplog.at_level(logging.INFO):
            rcode = computeMAICOS(system, logging.getLogger("test_logger"), ffonly=False)
        line = next(line for line in caplog.text.splitlines() if "Files to be computed:" in line)
        assert "TotalDensity" in line
        assert "DiporderWater" in line
        assert rcode == RCODE_ERROR
