# Using pulsar timing software to get constraints on P and Pdot

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
