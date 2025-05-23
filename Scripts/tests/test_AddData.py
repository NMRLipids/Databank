"""
`test_adddata.py` perform regression testing of adding-data functionality

NOTE: globally import of DatabankLib is **STRICTLY FORBIDDEN** because it 
      breaks the substitution of global path folders
"""

import os
import subprocess
from tempfile import TemporaryDirectory
import pytest


# run only without mocking data
pytestmark = pytest.mark.nodata

@pytest.fixture(scope="module")
def tmpWorkDir():
    with TemporaryDirectory(prefix="dbtestWD_", dir=os.path.dirname(__file__)) as wdir:
        print(f"Will use following directory for loadings: {wdir}")
        yield wdir


@pytest.fixture(scope="module")
def tmpOutDir():
    with TemporaryDirectory(prefix="Simulations.",
                            dir=os.path.join(os.path.dirname(__file__), "Data")
                            ) as wdir:
        print(f"Will use following directory for writting: {wdir}")
        yield wdir


class TestAddData:

    def _init(self):
        import DatabankLib
        if os.path.isfile(os.path.join(DatabankLib.NMLDB_DATA_PATH, '.notest')):
            pytest.exit("Test are corrupted. I see '.notest' file in the data folder.")
        self.exe = os.path.join(
            DatabankLib.NMLDB_ROOT_PATH, "Scripts",
            "BuildDatabank", "AddData.py")

    """
    Testing `AddData.py -h` behavior
    """

    def test_add_data_h(self):
        self._init()
        result = subprocess.run([
            self.exe,
            "-h"
        ], capture_output=True, text=True)
        print("STDOUT", result.stdout)
        print("STDERR:", result.stderr)
        assert result.returncode == 0

    """
    Testing `AddData.py -f <filename> -w <dirname> -o <dirname>` which should
    end correctly
    """

    @pytest.mark.parametrize(
            "infofn, debug", [("info566.yaml", False),
                              ("info566.yaml", True)])
    def test_add_data_addgood(self, infofn, debug, tmpWorkDir, tmpOutDir):
        self._init()
        fn = os.path.join(os.path.dirname(__file__), "Data", "info", infofn)
        runList = [
            self.exe,
            "-f",
            fn,
            "-w",
            tmpWorkDir,
            "-o",
            tmpOutDir
        ]
        if debug:
            runList.append("-d")
        result = subprocess.run(runList, capture_output=True, text=True)
        print("STDOUT", result.stdout)
        print("STDERR:", result.stderr)
        assert result.returncode == 0

    @pytest.mark.parametrize("infofn", ["info566_uf.yaml"])
    def test_add_data_fail(self, infofn, tmpWorkDir, tmpOutDir):
        self._init()
        fn = os.path.join(os.path.dirname(__file__), "Data", "info", infofn)
        result = subprocess.run([
            self.exe,
            "-f",
            fn,
            "-w",
            tmpWorkDir,
            "-o",
            tmpOutDir,
            "-d"
        ], capture_output=True, text=True)
        print("STDOUT", result.stdout)
        print("STDERR:", result.stderr)
        assert result.returncode == 1
