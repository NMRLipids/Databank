import os
import pytest
import logging

# Pytest HOOKS

def pytest_addoption(parser):
    parser.addoption(
        "--cmdopt", action="store", default="sim1", help="Two test groups: sim1|sim2|nodata"
    )

def pytest_collection_modifyitems(config, items):
    cmdopt = config.getoption("--cmdopt")
    run_only_sim = pytest.mark.skip(reason=f"Skipped by --cmdopt {cmdopt}")
    for item in items:
        if cmdopt not in item.keywords:
            item.add_marker(run_only_sim)

# Pytest GLOBAL FIXTURES

@pytest.fixture(autouse=True, scope="session")
def header_module_scope(request):
    cmdopt = request.config.getoption("--cmdopt")
    if cmdopt == "sim1":
        sim = "Simulations.1"
    elif cmdopt == "sim2":
        sim = "Simulations.2"
    elif cmdopt == "nodata":
        sim = None
    else:
        pytest.exit(f"Unknown --cmdopt {cmdopt}")
    os.environ["NMLDB_DATA_PATH"] = os.path.join(os.path.dirname(__file__), "Data")
    if sim is not None:
        os.environ["NMLDB_SIMU_PATH"] = os.path.join(os.path.dirname(__file__), "Data", sim)
    import DatabankLib
    print("DBG: Mocking Data path: ", DatabankLib.NMLDB_DATA_PATH)
    print("DBG: Mocking Simulations path: ", DatabankLib.NMLDB_SIMU_PATH)
    yield
    print("DBG: Mocking completed")

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