#!/usr/bin/env python

import astropy.units as u
import numpy as np
import copy
import sys
import pint.config
import pint.gridutils
import pint.models.parameter as param
import pint.residuals
import pint.logging
from pint.fitter import DownhillWLSFitter
from pint.models.model_builder import get_model, get_model_and_toas
from pint.toa import get_TOAs
pint.logging.setup(level="INFO")

import matplotlib.pyplot as plt
import matplotlib
from astropy.visualization import quantity_support

quantity_support()
import scipy.stats

def main():
    # Load in a basic dataset
    if len(sys.argv[1:]) != 2:
        print("Must provide the ephemeris (.par) file and the TOA (.tim) file")
        sys.exit(1)
    else:
        parfile = sys.argv[1]
        timfile = sys.argv[2]

    m, t = get_model_and_toas(parfile, timfile)

    f = DownhillWLSFitter(t, m)
    # Find the best-fit (or just compute statistics if no fitter parameters)
    f.fit_toas()
    bestfit = f.resids.chi2
    print(f"Initial fit chi2: {bestfit}")
    print(f"Free parameters: {f.model.free_params}")
    print("Fitter summary:")
    print(f.get_summary(nodmx=True))

    # Set up the gridding
    nsig_grid_f0 = 5  # Number of sigma to search around nominal value
    npts_f0 = 150  # Number of grid points for this parameter
    nsig_grid_f1 = 3
    npts_f1 = 300

    print(f"Setting up F0 grid ({npts_f0} points around nominal value +\- {nsig_grid_f0}-sigma)")
    print(f"  F0 = {f.model.F0.quantity} +/- {nsig_grid_f0 * f.model.F0.uncertainty}")
    F0 = np.linspace(
        f.model.F0.quantity - nsig_grid_f0 * f.model.F0.uncertainty,
        f.model.F0.quantity + nsig_grid_f0 * f.model.F0.uncertainty,
        npts_f0,
    )

    print(f"Setting up F1 grid ({npts_f1} points around nominal value +\- {nsig_grid_f1}-sigma)")
    print(f"  F1 = {f.model.F1.quantity} +/- {nsig_grid_f1 * f.model.F1.uncertainty}")
    F1 = np.linspace(
        f.model.F1.quantity - nsig_grid_f1 * f.model.F1.uncertainty,
        f.model.F1.quantity + nsig_grid_f1 * f.model.F1.uncertainty,
        npts_f1,
    )

    # Change PEPOCH to be in the middle of the TOA time span.
    # This is standard practice and can often help reduce degeneracies.
    print(f"Changing to PEPOCH {f.toas.get_mjds().mean()}")
    f.model.change_pepoch(f.toas.get_mjds().mean())

    # Compute the grid, where all parameters are frozen and sampled at
    # grid points, and the fitter is used to compute the fit statistic.
    print("Starting grid search...")
    chi2grid = pint.gridutils.grid_chisq(f, ("F0", "F1"), (F0, F1))[0]

    # Dump the output to a CSV file
    grid_f0, grid_f1 = np.meshgrid(F0.value, F1.value)
    with open("chi2_grid.csv", "w") as chi2file:
        chi2file.write("F0,F1,chi2,delta_chi2\n")
        for en, i in enumerate(np.ndindex(chi2grid.shape)):
            chi2file.write(f"{grid_f0[i]},{grid_f1[i]},{chi2grid[i]},{chi2grid[i] - bestfit}\n")

    # Get some of the grid statistics and nominal best values
    delta_chi2_min_idx = np.argmin(chi2grid - bestfit)
    nominal_best_f0 = grid_f0.flatten()[delta_chi2_min_idx]
    nominal_best_f1 = grid_f1.flatten()[delta_chi2_min_idx]
    print(f"Nominal best F0: {nominal_best_f0} Hz")
    print(f"Nominal best P0: {1/nominal_best_f0} s")
    print(f"Nominal best F1: {nominal_best_f1} Hz/s")
    print(f"Nominal best P1: {-nominal_best_f1/nominal_best_f0**2} s/s")
    print(f"Chi2 @ min. delta-chi2: {chi2grid.flatten()[delta_chi2_min_idx]}")
    print(f"Min. delta-chi2: {(chi2grid.flatten()-bestfit)[delta_chi2_min_idx]}")
    print(f"Min Chi2: {chi2grid.min()}")

    # 1, 2, and 3 sigma confidence limits
    nsigma = np.arange(1, 4)
    # These are the CDFs going from -infinity to nsigma.  So subtract away 0.5 and double for the 2-sided values
    CIs = (scipy.stats.norm().cdf(nsigma) - 0.5) * 2
    print(f"Confidence intervals for {nsigma} sigma: {CIs}")
    # chi^2 random variable for 2 parameters
    rv = scipy.stats.chi2(2)
    # The ppf = Percent point function is the inverse of the CDF
    contour_levels = rv.ppf(CIs)
    print(f"Contour levels for {nsigma} sigma and 2 parameters: {contour_levels}")

    # Do the same for a 1 parameter case
    CIs = (scipy.stats.norm().cdf(nsigma) - 0.5) * 2
    print(f"Confidence intervals for {nsigma} sigma: {CIs}")
    # chi^2 random variable for 1 parameters
    rv = scipy.stats.chi2(1)
    contour_levels_1param = rv.ppf(CIs)
    print(f"Contour levels for {nsigma} sigma and 1 parameter: {contour_levels_1param}")


    plots = True
    if plots:
        # Plot the grid/contour results
        fig, ax = plt.subplots(figsize=(16, 9))
        # Just plot the values offset from the best-fit values
        twod = ax.contour(
            F0 - f.model.F0.quantity,
            F1 - f.model.F1.quantity,
            chi2grid - bestfit,
            levels=contour_levels,
            colors="b",
        )
        oned = ax.contour(
            F0 - f.model.F0.quantity,
            F1 - f.model.F1.quantity,
            chi2grid - bestfit,
            levels=contour_levels_1param,
            colors="g",
            linestyles="--",
        )
        ax.errorbar(
            0,
            0,
            xerr=f.model.F0.uncertainty.value,
            yerr=f.model.F1.uncertainty.value,
            fmt="ro",
        )
        ax.set_xlabel("$\Delta F_0$ (Hz)", fontsize=24)
        ax.set_ylabel("$\Delta F_1$ (Hz/s)", fontsize=24)
        ax.set_title(
                f"Origin: $F_0 =$ {f.model.F0.quantity}   $F_1 =$ {f.model.F1.quantity}",
                fontsize=16
        )
        twod_artists, _ = twod.legend_elements()
        oned_artists, _ = oned.legend_elements()
        ax.legend(
            [twod_artists[0], oned_artists[0]],
            ["Joint 2D fit", "Single-parameter Fit"],
            fontsize=18, loc="upper left"
        )
        plt.savefig("gpm_f0f1_contours.png", format="png", bbox_inches="tight")
        plt.clf()


        # Now plot the actual chi^2 grid, with the contours overlaid
        delta_f0 = (F0[1] - F0[0]).value
        delta_f1 = (F1[1] - F1[0]).value
        f0_scale = 1e-9
        f0_scale_lbl = r"$\times 10^{{-9}}$"
        f1_scale = 1e-17
        f1_scale_lbl = r"$\times 10^{{-17}}$"


        im = plt.imshow(
            chi2grid - bestfit,
            origin="lower",
            interpolation="bicubic",
            cmap=plt.get_cmap("cubehelix_r"),
            norm=matplotlib.colors.SymLogNorm(10, vmax=1000),
            aspect="auto",
            extent=(
                    (F0[0].value - delta_f0 / 2 - f.model.F0.value) / f0_scale,
                    (F0[-1].value + delta_f0 / 2 - f.model.F0.value) / f0_scale,
                    (F1[0].value - delta_f1 / 2) / f1_scale,
                    (F1[-1].value + delta_f1 / 2) / f1_scale,
                ),
        )
        cs = plt.contour(
            (F0 - f.model.F0.quantity) / f0_scale,
            (F1 - f.model.F1.quantity) / f1_scale,
            chi2grid - bestfit,
            levels=contour_levels,
            colors=["k", "b", "r"],
            linewidths=1,
        )
        plt.grid(ls=":")
        plt.axvline(0, ls="--", lw=0.7, color="0.8")
        plt.axhline(0, ls="--", lw=0.7, color="0.8")
        plt.ylim(-2.1, 2.1)
        plt.xlim(-2.1, 2.1)
        plt.xlabel(
                rf"$F_0 - F_0^{{\rm init}}$ ({f0_scale_lbl} {f.model.F0.units})"
        )
        plt.ylabel(
                rf"$F_1$ ({f1_scale_lbl} {f.model.F1.units})"
        )
        cbar = plt.colorbar(im, label="Joint 2D fit $\Delta\chi^2$")
        cbar.add_lines(cs)
        for j, lab in enumerate(['$1\sigma$','$2\sigma$','$3\sigma$']):
            cbar.ax.text(-0.5, contour_levels[j], lab, ha='center', va='center')
        plt.savefig("gpm_f0f1_grid.png", format="png", bbox_inches="tight")


if __name__ == "__main__":
    main()
