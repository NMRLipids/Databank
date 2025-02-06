"""
`test_adddata.py` perform regression testing of adding-data functionality
"""

from unittest import mock
import os
import subprocess
from tempfile import TemporaryDirectory
import pytest
import DatabankLib


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


@pytest.fixture(autouse=True, scope="module")
def header_module_scope():
    _rp = os.path.join(os.path.dirname(__file__), "Data", "Simulations.1")
    with mock.patch.object(DatabankLib, "NMLDB_SIMU_PATH", _rp):
        print("DBG: Mocking simulation path: ", DatabankLib.NMLDB_SIMU_PATH)
        yield
    print("DBG: Mocking completed")


class TestAddData:

    exe = os.path.join(DatabankLib.NMLDB_ROOT_PATH, "Scripts",
                       "BuildDatabank", "AddData.py")

    """
    Testing `AddData.py -h` behavior
    """

    def test_add_data_h(self):
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
