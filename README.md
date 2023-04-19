# GPMTransient
Repository of code and data supporting the newly discovered active long-period radio transient

## Position

- **RA**: 18:39:02.0 (18.6505 hrs, 279.7583 deg)
- **Dec**: -10:31:49.5 (-10.5304 deg)
- **Uncertainty**: 0.15"

## Analyses

In all of the below, if there are `make` commands, you can use the `-j` command to parallelise execution.

<details>
<summary><b>Prepare lightcurves <i>without</i> barycentring</b></summary>

```
cd lightcurves_nobary
make lightcurves
cd ..
```

Expected output: `*_lightcurve.txt` files

</details>

<details>
<summary><b>Combining same-pulse lightcurves (and lag analysis)</b></summary>

```
cd lag_analysis
python lag_analysis.py
cd ..
```

Expected output: `*_lightcurve_mod.txt` files

</details>

<details>
<summary><b>Get the non-barycentred TOAs</b></summary>

```
cd toas_nobary
make all_toas_mod.tim
cd ..
```

Expected output: `all_toas_mod.tim`

</details>

<details>
<summary><b>Grid search in F0 and F1</b></summary>

```
cd P_Pdot_diagram
python grid_search_f_fdot.py dofit.par ../toas_nobary/all_toas_mod.tim
cd ..
```

</details>
