"""
`test_api` contains tests of API those functions, which doesn't require making
MDAnalysis Universe and recomputing something. These functions just read README
files and precomputed JSON data.

Test data is stored in `./Data/Simulations.2`

-------------------------------------------------------------------------------
Currently, we do mocking of NMLDB_SIMU_PATH and it requires to be done once
during pytest session. That's why we name this file "test2_api.py" because
it avoids `pytest` from executing it together with other "test_*.py" files.
Probably the problem can be solved using importlib and updating import specs.

NOTE: globally import of DatabankLib is **STRICTLY FORBIDDEN** because it 
      breaks the substitution of global path folders
"""

import os
import sys
import pytest

# run only on sim2 mocking data
pytestmark = pytest.mark.sim2

# for vector comparisons with np.testing.assert_allclose
MAXRELERR_COMPARE_THRESHOLD = 1e-2

@pytest.fixture(scope="module")
def systems():
    import DatabankLib
    from DatabankLib.core import initialize_databank
    if os.path.isfile(os.path.join(DatabankLib.NMLDB_DATA_PATH, '.notest')):
        pytest.exit("Test are corrupted. I see '.notest' file in the data folder.")
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
                         [(281, [0.261, 0.415, 0.609]),
                          (566, [0.260, 0.405, 0.597]),
                          (243, [0.281, 0.423, 0.638]),
                          (86,  [0.264, 0.419, 0.623])])
def test_GetFormFactorMin(systems, systemid, result):
    from DatabankLib.databankLibrary import GetFormFactorMin
    import numpy as np
    sys0 = systems.loc(systemid)
    ffl = GetFormFactorMin(sys0)
    np.testing.assert_allclose(
            np.array(ffl[:3]),
            np.array(result),
            rtol=MAXRELERR_COMPARE_THRESHOLD,
            err_msg=(
                     "Problem in FFMIN comparison:\n"
                  + f"Computed: {str(ffl[:3])} \n"
                  + f"Pre-computed: {str(result)}"
                ),
    )


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
         (787, ['TOCL/P3', 'POPC/P', 'POPE/P'], ['M_G13P1_M', 'M_G3P2_M', 'M_G3P2_M']),
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
         (243, 'nonExisting/P8', "was not found in the system")])
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


@pytest.mark.xfail(reason="Completely deprecated function", run=True,
                   raises=NotImplementedError, strict=True)
@pytest.mark.parametrize("systemid, lipid, result",
                         [(281, ['POPC'], [134])])
def test_raise_loadMappingFile(systems, systemid, lipid, result):
    from DatabankLib.databankLibrary import loadMappingFile
    sys0 = systems.loc(systemid)
    i = 0
    for lip in lipid:
        mpf = loadMappingFile(sys0['COMPOSITION'][lip]['MAPPING'])
        assert len(mpf) == result[i]
        i += 1


@pytest.mark.parametrize(
        "systemid, lipid, result",
        [(281, ['POPC/P'], ['M_G3P2_M']),
         (566, ['POPC/P31', 'CHOL/C1'], ['M_G3P2_M', 'M_C1_M']),
         (787, ['TOCL/P3', 'POPC/P', 'POPE/P'], ['M_G13P1_M', 'M_G3P2_M', 'M_G3P2_M']),
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
         (243, 'nonExisting/M_G1_M', "was not found in the system!")])
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
                         [(281, 1, 4.19996),
                          (243, 1, 4.25947)])
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


@pytest.mark.parametrize(
    "systemid, result",
    [
        (243, 0.7212884475213442),
        (86, 1.5018596337724872),
        (566, 1.1740608659926115),
    ],
)
def test_GetEquilibrationTimes(systems, systemid, result):
    from DatabankLib.databankLibrary import GetEquilibrationTimes
    from DatabankLib.settings.molecules import lipids_set

    sys0 = systems.loc(systemid)
    eq_times = GetEquilibrationTimes(sys0)
    print("\n========\n", eq_times, "\n=========\n")
    lips = list(set(sys0["COMPOSITION"].keys()).intersection(lipids_set.names))
    for lip in lips:
        if lip in ["CHOL", "DCHOL"]:
            continue  # for them eq_times are not computed
        assert lip in eq_times.keys()
        assert result == eq_times[lip]
    
@pytest.mark.xfail(reason="EQtimes were not computed", run=True,
                   raises=FileNotFoundError, strict=True)
def test_GetEquilibrationTimes_fail(systems):
    sys0 = systems.loc(787)
    from DatabankLib.databankLibrary import GetEquilibrationTimes
    GetEquilibrationTimes(sys0)
