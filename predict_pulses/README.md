# Predicting future TOAs

The following recipe can be used to generate TOAs between "now" (e.g. MJD 60066) and "some time in the future" (e.g. 60310).
The TOAs are written out in the TEMPO2 format of a `.tim` file.
The final steps (`cat` and `pintk`) are just to make sure that the predicted TOAs are sensible, i.e. consistent with the known ("input") TOAs.
(If they are not, the `--manual_offset` option may need to be used.)

```
# Produce the TOAs
#python predict_toas.py dofit.par ../toas_nobary/all_toas_mod.tim 60066 60310  --refit --manual_offset -640 --outparfile out.par
python predict_toas.py dofit.par ../toas_nobary/all_toas_mod.tim 60066 60310  --refit --outparfile out.par
# ^^^ This writes out 'out.tim' by default (can be changed with --outtimfile option).

# Check these alongside the observed TOAs to make sure all is well
cat ../toas_nobary/all_toas_mod.tim out.tim > all.tim
pintk out.par all.tim
```

# Filtering the predicted TOAs against some other criteria

Produce a list of GPS times which meet the desired local hour angle criterion:

```
python filter_toas.py out.tim out.txt --HA_range 22.5 1.5
```
