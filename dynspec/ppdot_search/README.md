# Using pulsar timing software to get constraints on P and Pdot

**NB:** All `make` commands below can be safely (and efficiently) run with parallelism (the `-j` option).

We estimate TOAs by convolving the pulses with a Gaussian that is (approximately) as wide as the known pulse window.
The pulse TOAs will be the peak lag, and the width of this Gaussian will be the errors on the TOAs.

The lightcurves need to be

- dedispersed to infinite frequency (this happens by default in `dedisperse_dynspec.py`), and
- *not* barycentred.

To do this, fresh `yaml` files need to be made that explicitly turn barycentring off.
This is done via an `awk` command in [`Makefile`](Makefile), which is invoked with

```
make prepare_yaml_files
```

which puts the revised `yaml` files into this directory.
Now we are ready to make the lightcurves themselves, by running the `dedisperse_dynspec.py` script on the new `yaml` files, which can be invoked with

```
make prepare_lightcurves
```

Finally, the script for generating the TOAs can be run.
The script deals with one lightcurve (i.e. one pulse) at a time.
To generate all the TOAs, and compile them into a single file (`all_toas.tim`), run

```
make all_toas.tim
```

This command will also generate plots of the lightcurves along with their low-pass filtered equivalents and their TOAs, which can be visually inspected to make sure the results are sensible.

The `all_toas.tim` file is also added to this repo so that it can be version controlled.

### Setting the TOA error

The error assigned to the TOAs is set in the [`Makefile`](Makefile) with the variable name `TOA_ERR_SEC`.
If you would like to change this and rerun the TOA calculations, it is recommended to

1. change the value in the Makefile,
2. run `touch *_lightcurves.txt`,
3. run `make all_toas.tim`.
