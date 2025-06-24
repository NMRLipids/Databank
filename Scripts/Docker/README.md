## More information about NMRlipids Docker images 

Right now there are:

1. `nmrlipids/gromacs`: This is the base gromacs image, only contains gromacs.

Then we have two images extending this gromacs base:

2. `nmrlipids/github-runner`: this is meant to be used by standard github runners which are sensitive to user configuration, spefically the UID and GID of the user. Handling python environments is also done in a more standardized way which fixed some errors caused by github's runner configurations.
3. `nmrlipids/core`: features everything developers need to work on the project including Conda. 

