"""
Test module test_ipylib.py

Contains tests of the plotting subpackage for IPython/Jupyter notebooks.

-------------------------------------------------------------------------------
NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import importlib
import sys

import pytest

pytestmark = [pytest.mark.nodata, pytest.mark.min]


def test_ipylib_imports_without_ipython(monkeypatch):
    """IPython is not a dependency, so ipylib must be importable without it."""
    monkeypatch.setitem(sys.modules, "IPython", None)  # makes `import IPython` fail
    for name in [n for n in sys.modules if n.startswith("fairmd.lipids.ipylib")]:
        monkeypatch.delitem(sys.modules, name)

    ipylib = importlib.import_module("fairmd.lipids.ipylib")
    assert callable(ipylib.plot_simulation_FF)
    assert callable(ipylib.plot_simulation_OP)
