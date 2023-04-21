# Using pulsar timing software to get constraints on P and Pdot

```
awk -F '[ ]' '/barycentric/ { $$4 = "false" } /Input file/ {$5 = "../dynspec/" $5} {print}' < ../dynspec/1234567890.yaml > 1234567890.yaml
```

or, to run on all available `../dynspec/*.yaml` files,

```
make yaml_files
```

which puts the revised `yaml` files into this directory.
Now we are ready to make the lightcurves themselves, by running the `dedisperse_dynspec.py` script on the new `yaml` files, which can be invoked with

```
make lightcurves
```
