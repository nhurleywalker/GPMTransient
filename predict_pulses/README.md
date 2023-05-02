# Predicting future TOAs

The following recipe can be used to generate TOAs between "now" (e.g. MJD 60066) and "some time in the future" (e.g. 60310).
The TOAs are written out in the TEMPO2 format of a `.tim` file.
The final steps (`cat` and `pintk`) are just to make sure that the predicted TOAs are sensible, i.e. consistent with the known ("input") TOAs.

```
#python predict_toas.py dofit.par ../toas_nobary/all_toas_mod.tim 60066 60310  --refit --manual_offset -640 --outparfile out.par
python predict_toas.py dofit.par ../toas_nobary/all_toas_mod.tim 60066 60310  --refit --outparfile out.par
cat ../toas_nobary/all_toas_mod.tim out.tim > all.tim
pintk out.par all.tim
```

# Filtering the predicted TOAs against some other criteria


