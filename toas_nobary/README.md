# Using pulsar timing software to get constraints on P and Pdot

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