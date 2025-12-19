# 1D-FW-CADIS Neutron Transport
This repository contains the software and report satisfying requirements laid out by NPRE 555: Reactor Theory I for Computer Project III.
In this project, I attempted to implement the FW-CADIS method to k-eigenvalue problems using the Monte Carlo code I created in Computer Project I.

## Code Compiling
To build the code:
```
mkdir build && cd build/
cmake ..
make -j6
```
Tests were not created (but the should be...).

## Examples
Example files have been included for convenience.
From `examples`, simply navigate to one of the directories, call the built executable, then give the input file name as a command line argument.
Note that examples in `/analog` will *not* run without modification.

## Running Simulations
`1D-FW-CADIS` takes one command-line argument, which is the name of the input file.
Input files are `XML` files, and the structure of these can be seen in the included examples.
Note that one can include and aribtrary number of regions so long as they are:
1. sorted from smallest x-values to largest
2. non-overlapping

Error messages have been written should the user forget these requirements.

Additionally, the method for variance-reduction must be specified. 
When in double, try something a read the error message :-).

## Plotting Results
For convenience, a plotting `Python` script as been included.
To use, follow this template:
```
python3 plot_outfile.py {input_filename}
```
A convention has been created such that only the input filename needs to be passed to get results plotted, *not* the output `.csv` filename.
Additionally, if one runs simulations with `fw-cadis` as variance reduction, the `plot_sn.py` script will plot the results of the discrete ordinates sweeps (both forward and adjoint) and also plot the importance map.
All values are normalized to have a max value of 1.0.

## Report
The report must be created.
To use the `Makefile` supplied, navigate to `report`, then run `make`.

### Replicating Results
The files `submission_simulation.xml` and `challenge_prob.xml` were the only problems tested and reported on.
Run as described above.

