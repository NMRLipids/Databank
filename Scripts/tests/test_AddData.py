"""
`test_adddata.py` perform regression testing of adding-data functionality

NOTE: globally import of DatabankLib is **STRICTLY FORBIDDEN** because it 
      breaks the substitution of global path folders
"""

import os
import shutil
import subprocess
from tempfile import TemporaryDirectory
import pytest


# run only without mocking data
pytestmark = pytest.mark.adddata


@pytest.fixture(scope="module")
def tmp_work_dir():
    with TemporaryDirectory(prefix="dbtestWD_", dir=os.path.dirname(__file__)) as wdir:
        print(f"Will use following directory for loadings: {wdir}")
        yield wdir


class TestAddData:

    @classmethod
    def setup_class(cls):
        import DatabankLib
        if os.path.isfile(os.path.join(DatabankLib.NMLDB_DATA_PATH, '.notest')):
            pytest.exit("Test are corrupted. I see '.notest' file in the data folder.")
        cls.exe = os.path.join(
            DatabankLib.NMLDB_ROOT_PATH, "Scripts",
            "BuildDatabank", "AddData.py")
        cls.out_dir = DatabankLib.NMLDB_SIMU_PATH
        os.mkdir(cls.out_dir)

    @classmethod
    def teardown_class(cls):
        if os.path.exists(cls.out_dir):
            shutil.rmtree(cls.out_dir, ignore_errors=True)

    def test_add_data_h(self):
        """
        Testing `AddData.py -h` behavior
        """
        result = subprocess.run([
            self.exe,
            "-h"
        ], capture_output=True, text=True)
        print(result.stdout)
        assert "-o OUTPUT_DIR", "Expected -o OUTPUT_DIR in stdout"
        assert result.returncode == 0

    @pytest.mark.parametrize(
            "infofn, debug", [("info566.yaml", False),
                              ("info566.yaml", True)])
    def test_add_data_addgood(self, infofn, debug, tmp_work_dir, capsys):
        """
        Testing `AddData.py -f <filename> -w <dirname> -o <dirname>` which should
        end correctly
        """
        fn = os.path.join(os.path.dirname(__file__), "Data", "info", infofn)
        run_list = [
            self.exe,
            "-f",
            fn,
            "-w",
            tmp_work_dir,
            "-o",
            self.out_dir
        ]
        if debug:
            run_list.append("-d")
        result = subprocess.run(run_list, capture_output=True, text=True)
        print(result.stderr)
        if debug:
            assert '[DEBUG]' in result.stderr
        else:
            assert '[DEBUG]' not in result.stderr
        assert '[INFO]' in result.stderr
        self._check_new_readmy(capsys)
        assert result.returncode == 0

    @pytest.mark.parametrize("infofn", ["info566_uf.yaml"])
    def test_add_data_fail(self, infofn, tmp_work_dir, capsys):
        fn = os.path.join(os.path.dirname(__file__), "Data", "info", infofn)
        result = subprocess.run([
            self.exe,
            "-f",
            fn,
            "-w",
            tmp_work_dir,
            "-o",
            self.out_dir,
            "-d"
        ], capture_output=True, text=True)
        print(result.stderr)
        assert '[ERROR]' in result.stderr, "Expected an error message in stderr"
        assert result.returncode == 1

    def _check_new_readmy(self, capsys):
        from DatabankLib.core import initialize_databank
        ss = initialize_databank()
        captured = capsys.readouterr()
        if '!!README LOAD ERROR!!' in captured.err:
            pytest.fail("There was an error in the captured output: " + captured.err)
        for s in ss:
            print(s['DOI'])
        assert len(ss) > 0, "No systems found in the databank after adding data."
