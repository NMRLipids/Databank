"""
`test_api` contains tests of API those functions, which doesn't require making
MDAnalysis Universe and recomputing something. These functions just read README
files and precomputed JSON data.

Test data is stored in `./Data/Simulations.2`

-------------------------------------------------------------------------------
TODO think to repair the following:
Currently, we do mocking of NMLDB_SIMU_PATH and it requires to be done once
during pytest session. That's why we name this file "test2_api.py" because
it avoids `pytest` from executing it together with other "test_*.py" files.
Probably the problem can be solved using importlib and updating import specs.
"""

from unittest import mock
import os
import sys
import pytest
import logging


@pytest.fixture(scope="module")
def logger():
    logger = logging.getLogger("test_logger")
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)

    if not logger.handlers:  # Avoid adding multiple handlers during pytest runs
        logger.addHandler(handler)

    yield logger

    # TEARDOWN: clean up handlers after use
    for handler in logger.handlers:
        logger.removeHandler(handler)


@pytest.fixture(autouse=True, scope="module")
def header_module_scope():
    import DatabankLib
    _rp = os.path.join(os.path.dirname(__file__), "Data", "Simulations.2")
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
    assert len(systems) == 5


def test_print_README(systems, capsys):
    from DatabankLib.core import print_README
    sys0 = systems[0]
    print_README(sys0)
    output: str = capsys.readouterr().out.rstrip()
    fDOI = output.find('DOI:') != -1
    fTEMP = output.find('TEMPERATURE:') != -1
    assert fDOI and fTEMP


@pytest.mark.parametrize("systemid, result",
                         [(281, 64.722),
                          (566, 61.306),
                          (787, 78.237),
                          (243, 62.276),
                          (86,  60.460)])
def test_CalcAreaPerMolecule(systems, systemid, result):
    from DatabankLib.databankLibrary import CalcAreaPerMolecule
    sys0 = systems.loc(systemid)
    apm = CalcAreaPerMolecule(sys0)
    assert abs(apm - result) < 6e-4


@pytest.mark.parametrize("systemid, result",
                         [(281, 4142.234),
                          (566, 3923.568),
                          (787, 4694.191),
                          (243, 2241.920),
                          (86,  3869.417)])
def test_calcArea(systems, systemid, result):
    from DatabankLib.databankLibrary import calcArea
    sys0 = systems.loc(systemid)
    area = calcArea(sys0)
    assert abs(area - result) < 6e-4


@pytest.mark.parametrize("systemid, result",
                         [(281, 128),
                          (566, 128),
                          (787, 120),
                          (243, 72),
                          (86,  128)])
def test_GetNLipids(systems, systemid, result):
    from DatabankLib.databankLibrary import GetNlipids
    sys0 = systems.loc(systemid)
    nlip = GetNlipids(sys0)
    assert nlip == result


@pytest.mark.parametrize("systemid, result",
                         [(281, [0.264, 0.415, 0.614]),
                          (566, [0.263, 0.404, 0.604]),
                          (243, [0.281, 0.423, 0.638]),
                          (86,  [0.266, 0.421, 0.63])])
def test_GetFormFactorMin(systems, systemid, result):
    from DatabankLib.databankLibrary import GetFormFactorMin
    import numpy as np
    sys0 = systems.loc(systemid)
    ffl = GetFormFactorMin(sys0)
    err = ((np.array(ffl[:3]) - np.array(result))**2).sum()
    assert err < 1e-9


@pytest.mark.parametrize("systemid, result",
                         [(281, 31.5625),
                          (566, 31.0),
                          (787, 75.0),
                          (243, 39.7778),
                          (86,  27.75)])
def test_getHydrationLevel(systems, systemid, result):
    from DatabankLib.databankLibrary import getHydrationLevel
    sys0 = systems.loc(systemid)
    hl = getHydrationLevel(sys0)
    assert abs(hl - result) < 1e-4


@pytest.mark.parametrize("systemid, lipid, result",
                         [(281, ['POPC'], [1]),
                          (566, ['POPC', 'CHOL'], [.9375, .0625]),
                          (787, ['TOCL', 'POPC', 'POPE'], [0.25, 0.5, 0.25]),
                          (243, ['DPPC'], [1]),
                          (86,  ['POPE'], [1])])
def test_calcLipidFraction(systems, systemid, lipid, result):
    from DatabankLib.databankLibrary import calcLipidFraction
    sys0 = systems.loc(systemid)
    assert calcLipidFraction(sys0, 'SOPC') == 0  # absent lipid
    err = 0
    i = 0
    for lip in lipid:
        err += (calcLipidFraction(sys0, lip) - result[i])**2
        i += 1
    assert err < 1e-4


@pytest.mark.parametrize("systemid, result",
                         [(281, [-0.1610, -0.1217]),
                          (566, [-0.1714, -0.1142]),
                          (243, [-0.1764, -0.1784]),
                          (86,  [-0.1933, -0.1568])])
def test_averageOrderParameters(systems, systemid, result):
    from DatabankLib.databankLibrary import averageOrderParameters
    sys0 = systems.loc(systemid)
    sn1, sn2 = averageOrderParameters(sys0)
    assert (sn1-result[0])**2 + (sn2-result[1])**2 < 1e-5

# Tests behavior when averageOrderParameters cannot find calculated OP data


@pytest.mark.parametrize("systemid", [787], )
def test_raises_averageOrderParameters(systems, systemid):
    from DatabankLib.databankLibrary import averageOrderParameters
    sys0 = systems.loc(systemid)
    with pytest.raises(FileNotFoundError) as exc_info:
        # sn1, sn2 = ..
        averageOrderParameters(sys0)
    assert 'OrderParameters.json' in str(exc_info.value)


@pytest.mark.parametrize(
        "systemid, lipid, result",
        [(281, ['POPC/P'], ['M_G3P2_M']),
         (566, ['POPC/P31', 'CHOL/C1'], ['M_G3P2_M', 'M_C1_M']),
         (787, ['TOCL/P3', 'POPC/P', 'POPE/P'], ['M_G13P2_M', 'M_G3P2_M', 'M_G3P2_M']),
         (243, ['DPPC/P8'], ['M_G3P2_M']),
         (86,  ['POPE/P8'], ['M_G3P2_M'])])
def test_getUniversalAtomName(systems, systemid, lipid, result):
    from DatabankLib.databankLibrary import getUniversalAtomName
    sys0 = systems.loc(systemid)
    i = 0
    for lipat in lipid:
        lip, atom = tuple(lipat.split('/'))
        uname = getUniversalAtomName(sys0, atom, lip)
        assert uname == result[i]
        i += 1

# Test fail-behavior of getUniversalAtomName


@pytest.mark.parametrize(
        "systemid, lipat, result",
        [(243, 'DPPC/nonExisting', "Atom was not found"),
         (243, 'nonExisting/P8', "Mapping file was not found")])
def test_bad_getUniversalAtomName(systems, systemid, lipat, result, capsys):
    from DatabankLib.databankLibrary import getUniversalAtomName
    sys0 = systems.loc(systemid)
    lip, atom = tuple(lipat.split('/'))
    uname = getUniversalAtomName(sys0, atom, lip)
    output = capsys.readouterr().err.rstrip()
    assert result in output
    assert uname is None


@pytest.mark.parametrize("systemid, lipid, result",
                         [(243, 'DPPC', "44ea5"),
                          (787, 'TOCL', "78629")])
def test_getAtoms(systems, systemid, lipid, result):
    from DatabankLib.databankLibrary import getAtoms
    sys0 = systems.loc(systemid)
    atoms = getAtoms(sys0, lipid).split()
    atoms = ",".join(sorted(atoms))
    import hashlib
    md5_hash = hashlib.md5()
    md5_hash.update(atoms.encode('ascii'))
    hx = md5_hash.hexdigest()[:5]
    assert hx == result


@pytest.mark.parametrize("systemid, lipid, result",
                         [(281, ['POPC'], [134]),
                          (566, ['POPC', 'CHOL'], [134, 74]),
                          (787, ['TOCL', 'POPC', 'POPE'], [248, 134, 125]),
                          (243, ['DPPC'], [130]),
                          (86,  ['POPE'], [125])])
def test_loadMappingFile(systems, systemid, lipid, result):
    from DatabankLib.databankLibrary import loadMappingFile
    sys0 = systems.loc(systemid)
    i = 0
    for lip in lipid:
        mpf = loadMappingFile(sys0['COMPOSITION'][lip]['MAPPING'])
        assert len(mpf) == result[i]
        i += 1


@pytest.mark.xfail(reason="Improper file name", run=True,
                   raises=FileNotFoundError, strict=True)
def test_raise_loadMappingFile():
    from DatabankLib.databankLibrary import loadMappingFile
    mpf = loadMappingFile('file-doesnot-exist')
    print(mpf)


@pytest.mark.parametrize(
        "systemid, lipid, result",
        [(281, ['POPC/P'], ['M_G3P2_M']),
         (566, ['POPC/P31', 'CHOL/C1'], ['M_G3P2_M', 'M_C1_M']),
         (787, ['TOCL/P3', 'POPC/P', 'POPE/P'], ['M_G13P2_M', 'M_G3P2_M', 'M_G3P2_M']),
         (243, ['DPPC/P8'], ['M_G3P2_M']),
         (86,  ['POPE/P8'], ['M_G3P2_M'])])
def test_simulation2universal_atomnames(systems, systemid, lipid, result):
    from DatabankLib.databankLibrary import simulation2universal_atomnames
    sys0 = systems.loc(systemid)
    i = 0
    for lipat in lipid:
        lip, atom = tuple(lipat.split('/'))
        sname = simulation2universal_atomnames(sys0, lip, result[i])
        assert sname == atom
        i += 1


@pytest.mark.parametrize(
        "systemid, lipat, result",
        [(243, 'DPPC/nonExisting', "was not found from mappingDPPCberger.yaml"),
         (243, 'nonExisting/M_G1_M', "Mapping file was not found")])
def test_bad_simulation2universal_atomnames(systems, systemid, lipat, result, capsys):
    from DatabankLib.databankLibrary import simulation2universal_atomnames
    sys0 = systems.loc(systemid)
    lip, atom = tuple(lipat.split('/'))
    sname = simulation2universal_atomnames(sys0, lip, atom)
    output = capsys.readouterr().err.rstrip()
    assert result in output
    assert sname is None


@pytest.mark.parametrize(
        "systemid, result",
        [(281, "resname POPC"),
         (566, "resname CHL or resname OL or resname PA or resname PC"),
         (787, "resname POPC or resname POPE or resname TOCL2"),
         (243, "resname DPPC"),
         (86,  "resname POPE")])
def test_getLipids(systems, systemid, result):
    from DatabankLib.databankLibrary import getLipids
    sys0 = systems.loc(systemid)
    gl = getLipids(sys0)
    assert gl == result


# TEST thickness calculation here because it is not trajectory-based, but JSON-based
"""
@pytest.mark.parametrize("systemid, result",
                         [(566, 0),
                          (787, 2),
                          (86,  0)])
def test_analyze_th(systems, systemid, result):
    import DatabankLib
    from DatabankLib.analyze import computeTH
    sys0 = systems.loc(systemid)
    rc = computeTH(sys0)
    assert rc == result
    if rc == DatabankLib.RCODE_ERROR:
        fn = os.path.join(
            DatabankLib.NMLDB_SIMU_PATH,
            sys0['path'], 'thickness.json')
        assert not os.path.isfile(fn)  # file is not created
"""


@pytest.fixture(scope='function')
def wipeth(systems):
    import DatabankLib
    # TD-FIXTURE FOR REMOVING THICKNESS AFTER TEST CALCULATIONS
    yield
    # TEARDOWN
    for sid in [243, 281]:
        sys0 = systems.loc(sid)
        fn = os.path.join(DatabankLib.NMLDB_SIMU_PATH, sys0['path'], 'thickness.json')
        try:
            os.remove(fn)
        except FileNotFoundError:
            pass
        except Exception as e:
            print(f"An error occured during teardown: {e}", file=sys.stderr)
            raise


@pytest.mark.parametrize("systemid, result, thickres",
                         [(281, 1, 4.18335),
                          (243, 1, 4.27262)])
def test_analyze_th(systems, systemid, result, wipeth, thickres, logger):
    import DatabankLib
    from DatabankLib.analyze import computeTH
    sys0 = systems.loc(systemid)
    rc = computeTH(sys0, logger)
    assert rc == result
    fn = os.path.join(
        DatabankLib.NMLDB_SIMU_PATH,
        sys0['path'], 'thickness.json')
    assert os.path.isfile(fn)
    with open(fn, 'r') as file:
        data = float(file.read().rstrip())
    assert abs(data - thickres) < 1e-5


@pytest.mark.parametrize("systemid, result",
                         [(281, None),
                          (243, None),
                          (566, 4.2576),
                          (787, None),
                          (86,  4.1327)])
def test_GetThickness(systems, systemid, result):
    from DatabankLib.databankLibrary import GetThickness
    sys0 = systems.loc(systemid)
    th = GetThickness(sys0)
    assert (th is None and result is None) or abs(float(th) - float(result)) < 1e-4


@pytest.mark.parametrize("systemid, result",
                         [(243, "0.7212"),
                          (86, "1.5018"),
                          (566, "1.174")])
def test_ShowEquilibrationTimes(systems, capsys, systemid, result):
    from DatabankLib.databankLibrary import ShowEquilibrationTimes
    from DatabankLib.settings.molecules import lipids_dict
    sys0 = systems.loc(systemid)
    _ = ShowEquilibrationTimes(sys0)
    captured = capsys.readouterr()
    print("\n========\n", captured, "\n=========\n")
    lips = list(set(sys0['COMPOSITION'].keys()).intersection(lipids_dict.keys()))
    for lip in lips:
        if lip in ["CHOL", "DCHOL"]:
            continue  # for them eq_times are not computed
        assert lip in captured.out
    assert result in captured.out


@pytest.mark.xfail(reason="EQtimes were not computed", run=True,
                   raises=FileNotFoundError, strict=True)
def test_ShowEquilibrationTimes_fail(systems):
    sys0 = systems.loc(787)
    from DatabankLib.databankLibrary import ShowEquilibrationTimes
    ShowEquilibrationTimes(sys0)
