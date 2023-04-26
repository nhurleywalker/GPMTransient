# Get TOAs for each barycentred pulse

```
make all_toas.tim
```

This command will also generate plots of the lightcurves along with their low-pass filtered equivalents and their TOAs, which can be visually inspected to make sure the results are sensible.

### Setting the TOA error

The error assigned to the TOAs is set in the [`Makefile`](Makefile) with the variable name `TOA_ERR_SEC`.
If you would like to change this and rerun the TOA calculations, it is recommended to

1. change the value in the Makefile,
2. run `touch *_lightcurve.txt`,
3. run `make all_toas.tim`.
