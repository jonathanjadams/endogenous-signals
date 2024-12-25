# endogenous-signals

This repository contains the replication package for the following paper:

* Jonathan J Adams "Macroeconomic Models with Incomplete Information and Endogenous Signals", (2024)

# Required Software and Packages

Computation was carried out in `Matlab`.  All results are confirmed on `Windows` using `Matlab` version `2021a`

## Replication Instructions

1. Ensure the `SIGNAL_OP Package` subfolder is included in your `Matlab` path
2. By default, plots will be generated but not saved.  If you want to automatically save the generated plots, create a subfolder `graphs`, and edit the `MMIIES_replication_6` file to include `plot_saving = 1`.
3. Run `MMIIES_replication_6`


# Repository Contents

`Matlab` code:

* `MMIIES_replication_6`: This is the main replication file.  Use this file to reproduce the plots from the paper. All other files are dependencies or sub-dependencies.
* `MMIIES_beauty_contest`: This file solves the indeterminacy regions in the beauty contest model.
* `MMIIES_singleton`: This file compares solutions to the Singleton model
* `MMIIES_confoundingdynamics`: This file computes the indeterminacy region in the idiosyncratic confounding dynamics model, and the perturbation divergence in the baseline confounding dynamics model
* `MMIIES_singleton_ztran`: This file solves the Singleton model with the Han et al (2022) z-tran toolbox; it produces the output fils "ztran_singleton" and "ztran_singleton_time".  This file does not need to be run to reproduce the paper figures.  But to make comparisons about algorithm solution times using your own computer, you will need to run this program.
* `beauty_contest_discriminant`: sub-routine used in `MMIIES_beauty_contest`
* `cdi_discriminant`: sub-routine used in `MMIIES_confoundingdynamics` 

Folder:

* `SIGNAL_OP Package`: subroutines for solving models by Signal Operator Iteration

Matrices:

* `nimark_pirfs`: Impulse response functions from the Nimark algorithm (not publicly available)
* `ztran_singleton`: Output of impulse response functions from `MMIIES_singleton_ztran`
* `ztran_singleton_time`: Model solution time from `MMIIES_singleton_ztran`


