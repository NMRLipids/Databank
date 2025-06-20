# Code style and naming conventions

We use `flake8` as a linter-check with `pep8-naming` conventions. 
So, we use `variable_name`, `ClassName`, and `method_name`. For global variables, we use `UPPER_CASE_NAMES`.
Sometimes, this convention can be broken, but a rule should be written in the end to allow vscode to ignore the breach.

We use type hints and `typing` library for hints like `List[int]`.

We inherit `collections.abc` classes if we need to construct a collection with an extended behavior.

We are using automatic testing with `pytest` please consult (test-README)[Scripts/tests/README.md]. We use both unit and regression testing, however the coverage of code with unit-tests remains relatively low.

# Repository rules

1. Please create an issue, and then pull-request, not vice versa.
2. Please keep main branch in your fork updated so that it can be rebased (we prefer rebasing over merging).
3. Once you address an issue, please refer it in the commit message by using a phrase such a 'Partially fixes #000'.
4. You must make sure that the tests don't break before you make a PR.
5. We don't oblige contributors to squash commits. We appreciate if the commit has a considerable size.
6. We require at least one review from organisation member to accept a pull-request.
7. Check representation of yourself in `AUTHORS.md` file.
8. You are always welcome to participate in repository discussions and NMRlipids community events to develop code together.


# PATH handling in our code

As the code operates with an over-filesystem database, universal path handling is crucial. To construct paths, we use global variables of the type NMLDB_**XXXX**_PATH, where **XXXX** could be:
- **ROOT** points by default to the root folder of the repository. Can be set-up and *should not* be used for Data-related paths. Once set alone, all nested paths are recreated from it.
- DATA points by default to `./Data` folder. Can be used to point to somewhere inside `Data` but not for molecules, simulations, and experiments. Currently is used for Rankings as well.
- MOL points by default to `./Data/Molecules` folder, from where molecule lists are initialized. Is used to get access to something in molecule-folders.
- SIMU points by default to `./Data/Simulations` folder. Is used to get access to a certain simulation.
- EXP points by default to `./Data/Experiments` folder. Is used to get access to a certain experiment.

We currently construct paths using `os.path.join(a,b)`.

# Documentation generation

Our documentation is automatically generated and deployed to https://nmrlipids.github.io. 
We use `sphinx` with Read-The-Docs plugin to fullfill documentation pages. 
Please refer [this page about Sphinx RTD](https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html) to guide docstrings styles.


# Data handling

NMRlipids Databank separates codespace from [Data](https://github.com/NMRLipids/BilayerData) since June 2025 (v.1.1.0). Data contribution rules are moved there accordingly.