# B*B simulation of the beam-beam interaction at LHC during van der Meer scan luminosity calibration

## Simulation
The `B*B` (B-star-B or `BxB`) simulation splits the "kicked" bunch into
individual "macro-particles" with weights, and traces them in the
accelerator. Every macro-particle propagates between the interaction points
(IPs) according to the specified betatron phase advances. At every IP it is
influenced by the electromagnetic interaction with the opposite "kicker"
bunch. The original shape of each bunch can be specified as an arbitrary
weighted sum of Gaussians with common mean, round or elliptical in X-Y
projection. The number of IPs is also arbitrary. In the end of the simulation
the modifications of the bunch overlap integrals (ie. of the luminosities) due
to beam-beam are reported.

## Usage
The main function is called `beam_beam`. There are 4 ways to call it.

### From R language
See the corresponding R package
https://github.com/balagura/BxB which contains the same code callable from R.

### From python 
One needs to compile a shared library `libBxB.so` by 
```
make libBxB.so
```
The example of the python program with explanations is given in
`BxB_example.py`. `cpython` is used to call C from python.

### From C++ program
Please, include "bb.hh" header (which in turn includes
bb_C_CPP.h with common C and C++ declarations). It contains definitions of
`Kicked`, `Kickers` and `Sim` structures for the kicked bunch, kickers and the
simulation parameters, respectively. They should be filled before calling the
main function `beam_beam(kicked, kickers, sim, &summary, quiet)`. The output
is a matrix `vector<vector<BB_Summary_Per_Step_IP> > summary`, with
`summary[ip][step]` containing the simulation results. The library `libBxB.so`
created by `make libBxB.so` and containing `beam_beam()` should be included (as
`-lBxB`) when making a C++ executable.

### Using configuration text file
The standalone C++ program `beam_beam` from `beam_beam.cpp` runs the simulation with 
the input parameters from the `Config` (`config.hh/cpp`) file given as a first argument. 
To build: 
```
make beam_beam
```
One example of the configuration is given in `BxB_example.txt` with the explanation of
all parameters and the simulation in general.

In all 4 cases the same C++ function `beam_beam()` from `bb.cpp` is called
internally. The input and output data corresponding to the `Kicked`,
`Kickers`, `Sim` and the matrix of `BB_Summary_Per_Step_IP` structures,
respectively, are also the same.

In addition to `BB_Summary_Per_Step_IP` matrix, the simulation can print more
detailed information **in separate files in the `output_dir`** subdirectory
specified in `Sim`. What is printed is controlled by another parameter **`output`**
which is a white-space separated list of available **options**.

## Installation and documentation
To start working with the simulation, one needs 
```
git clone https://github.com/balagura/beam-beam-simulation-for-vdM-scans-at-LHC/
make <desired target>
```
and run the program for the appropriately modified set of
input parameters. Their full description can be found in `BxB_example.py`,
`BxB_example.txt` or in `bb.hh`. Alternatively, one can install R package from
https://github.com/balagura/BxB and type in R: 
```
library(BxB)
?beam_beam
?kicked
?kickers
?sim
```

## Old C++ program for validation of `B*B` model with single-IP single-Gaussian round bunches
It was used to validate `B*B` with various kick models and to compare it 
  with the simulation used at LHC before, in 2012-2018.
  It can be build with 
  ```
  make beam_beam_model_validation
  ```
The minimal (possibly, outdated) documentation prepared at that time can be obtained as a PDF
document by typing `make` in the subdirectory `doc`.
