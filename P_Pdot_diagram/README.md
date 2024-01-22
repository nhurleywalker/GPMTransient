# P-Pdot analysis and plots

## P-Pdot diagram

```
python plot_P_Pdot.py
```

## P-Pdot search

First, generate a chi-squared grid on P and Pdot.
```
python grid_search_f_fdot.py
```
This requires the Python `pint-pulsar` module, and may take a while to run.
This outputs a file called `chi2_grid.csv`.
If you would like to skip this step, you can copy the pre-generated `chi2_grid.csv` file from the `upload_materials` folder into this folder.

## Generate a plot of the chi-squared as a function of P and Pdot

```
python plot_P_Pdot_search.py
```
