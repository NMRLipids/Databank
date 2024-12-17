*TODO: move to CI/CD, remove build from the repo!*

To recreate documentation, install:

```
apt install python3-sphinx
pip install myst-parser
pip install sphinx-rtd-theme
```

**Note!** You should be already in the environment suitable for working with NMRLipids Databank project.

Then run:
```
make html
```