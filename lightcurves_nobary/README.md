# Using pulsar timing software to get constraints on P and Pdot

**NB:** All `make` commands below can be safely (and efficiently) run with parallelism (the `-j` option).

We estimate TOAs by convolving the pulses with a Gaussian that is (approximately) as wide as the known pulse window.
The pulse TOAs will be the peak lag, and the width of this Gaussian will be the errors on the TOAs.

The lightcurves need to be

- dedispersed to infinite frequency (this happens by default in `dedisperse_dynspec.py`), and
- *not* barycentred.

To do this, fresh `yaml` files need to be made that explicitly turn barycentring off.
The yaml files also need to include a change to the relative path to the dynamic spectrum data.
Both of these can done via an `awk` command, e.g.

```
awk -F '[ ]' '/barycentric/ { $$4 = "false" } /Input file/ {$5 = "../dynspec/" $5} {print}' < ../dynspec/1234567890.yaml > 1234567890.yaml
```

or, to run on all available `../dynspec/*.yaml` files,

```
make prepare_yaml_files
```

which puts the revised `yaml` files into this directory.
Now we are ready to make the lightcurves themselves, by running the `dedisperse_dynspec.py` script on the new `yaml` files, which can be invoked with

```
make prepare_lightcurves
```
