# MCSlab
This repository contains the software and report satisfying requirements laid out by NPRE 555: Reactor Theory I for Computer Project I.

## Code Compiling
To build the code, navigate to `MCSlab/`.
Then:
```
mkdir build && cd build/
cmake ..
make
```
To ensure everything is running correctly, also run the tests.
Within `build/`, run:
```
./run_tests
```

## Examples
Example files have been included for convenience.
To run all, from the root directory:
```
cd MCSlab/examples
```
then:
```
bash run_examples.sh
```

## Running Simulations
`MCSlab` takes one command-line argument, which is the name of the input file.
Input files are `XML` files, and the structure of these can be seen in the included examples.
Note that one can include and aribtrary number of regions so long as they are:
1. sorted from smallest x-values to largest
2. non-overlapping

Error messages have been written should the user forget these requirements.

## Plotting Results
For convenience, a plotting `Python` script as been included.
To use, follow this template:
```
python3 plot_simulations {input_filename}
```
A convention has been created within `MCSlab` such that only the input filename needs to be passed to get results plotted, *not* the output `.csv` filename.
The user is encouraged to view the bash script included in `examples/` for troubleshooting.

## Report
The report must be created.
To use the `Makefile` supplied, navigate to `report`, then run `make`.

### Replicating Results
The random number generator within `MCSlab` is randomly seeded for each generation, so perfect replication is not possible. 
However, users are able to run the same input files because all are included within `examples/`.
Runnning the provided bash script will reproduce the results.
