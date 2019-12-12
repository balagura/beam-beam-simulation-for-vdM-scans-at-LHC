## Old code for validation of `BxB` model from version v1.0.
This subdirectory contains the old C++ program `beam_beam_model_validation.cpp` equivalent to `beam_beam.C` in v1.0. 
It was used for a validation of `B*B` model with single-IP single-Gaussian round bunches and to compare it with 
the old MAD-X'12 simultion used in LHC in 2012-2018.

There are many kick models which can be simulated by `beam_beam_model_validation.cpp` using `compare_models.sh`
```
precise
precise.minus.average
average
precise.minus.one.third.of.average
one.third.of.average
precise.minus.at.center
at.center
quadrupole
average.and.quadrupole
```
with adiabatic (n_transitional_turns=1000) or abrupt (n_transitional_turns=1) beam-beam switch on.

An old prescription of the beam-beam corrections at LHC in 2012-2018 was essentially equivalent to `average.and.quadrupole` 
kick model, ie. the electric field was the sum of the constant dipole field (`average`) and the field proportional to
the displacement from the kicker center in X and Y (`quadrupole`).

`Precise`, `precise.minus.average` and `average` are the kick models which are only kept in `BxB` v2.0. 

In `precise.minus.one.third.of.average` and `one.third.of.average` only 1/3 of the average kick is taken. 
In `precise.minus.at.center` and `at.center` the field value at the kicked bunch center is taken as a 
constant dipole field instead of the average over the full bunch.
.

Executing
```
./compare_models.sh
```
will build `beam_beam_model_validation` (this can also be done by `make beam_beam_model_validation`),
simulate all kick models with adiabatic and abrupt switch on, will analyze the data and produce various plots in
`bb_MAD-X_2012_*/pdf` subdirectories using R code from `comparison.R`. The main plots are in 
`bb_MAD-X_2012_precise.minus.*/pdf/gg.cor.dipole.subtr.and.added.pdf`.

## Minimal documentation
Old (possibly, outdated) documentation prepared at that time can be obtained as a PDF
document by typing `make` in the subdirectory `doc`.
