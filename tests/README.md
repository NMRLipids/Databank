# Testing strategy

Test function requiring `Simulation.1` and `Simulation.2` path patches should be run 
independently because otherwise mocking DatabankLib module interfere between modules.

So, simply run

```
 pytest -vs
```

to run everything network-dependent from `./Data/Simulation.1` subfolder. And then: 

```
 pytest test2_api.py -vs
```

to test API function on pre-computed toy databank.

# Explanation

In `./Data/Simulation.1` we collected a toy databank with only README files, so that testing 
functions can download everything, perform analyzis over the trajectory and write files. This
behavior corresponds to the situation when _de novo_ analysis is run after adding a new system.

In `./Data/Simulation.2`, in contrary, we put completed Databank analysis files and we will test 
only those API part on this toy-Databank, which doesn't require MDAnalysis objects to create. Some
function requires only JSON files to exist - for them we prepared this folder.

TODO: complete