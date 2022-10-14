#from astropy.table import Table
#from astropy.coordinates import SkyCoord
#from astropy import unit as u
#import scipy.ndimage
#from scipy.ndimage.filters import gaussian_filter
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import markers
from matplotlib.path import Path
import matplotlib.font_manager
from matplotlib import rc

# Nature requires sans-serif fonts
plt.rcParams.update({
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

cm = 1/2.54  # centimeters in inches

# https://stackoverflow.com/questions/26686722/align-matplotlib-scatter-marker-left-and-or-right
# To align the limits with their real values
def align_marker(marker, halign='center', valign='middle',):
    """
    create markers with specified alignment.

    Parameters
    ----------

    marker : a valid marker specification.
      See mpl.markers

    halign : string, float {'left', 'center', 'right'}
      Specifies the horizontal alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'center',
      -1 is 'right', 1 is 'left').

    valign : string, float {'top', 'middle', 'bottom'}
      Specifies the vertical alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'middle',
      -1 is 'top', 1 is 'bottom').

    Returns
    -------

    marker_array : numpy.ndarray
      A Nx2 array that specifies the marker path relative to the
      plot target point at (0, 0).

    Notes
    -----
    The mark_array can be passed directly to ax.plot and ax.scatter, e.g.::

        ax.plot(1, 1, marker=align_marker('>', 'left'))

    """

    if isinstance(halign, str): #, unicode)):
        halign = {'right': -1.,
                  'middle': 0.,
                  'center': 0.,
                  'left': 1.,
                  }[halign]

    if isinstance(valign, str): #, unicode)):
        valign = {'top': -1.,
                  'middle': 0.,
                  'center': 0.,
                  'bottom': 1.,
                  }[valign]

    # Define the base marker
    bm = markers.MarkerStyle(marker)

    # Get the marker path and apply the marker transform to get the
    # actual marker vertices (they should all be in a unit-square
    # centered at (0, 0))
    m_arr = bm.get_path().transformed(bm.get_transform()).vertices

    # Shift the marker vertices for the specified alignment.
    m_arr[:, 0] += halign / 2
    m_arr[:, 1] += valign / 2

    return Path(m_arr, bm.get_path().codes)

#from labellines import *

# Define function for string formatting of scientific notation
# Adapted from:
# https://stackoverflow.com/questions/18311909/how-do-i-annotate-with-power-of-ten-formatting
def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if exponent is None:
        exponent = int(np.floor(np.log10(np.abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits

#    return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)
    return r"$10^{{{0:d}}}$".format(exponent)


def Pdot_from_P_B(P, B):
    ''' feed an array of periods, and a single value for magnetic field.
        returns an array of period derivatives the same length as the period array'''
    return ( (B / 6.4e19 )**2 ) / P

def Pdot_from_P_E(P, E):
    ''' feed an array of periods, and a single value for spin down luminosity.
        returns an array of period derivatives the same length as the period array'''
    return E * P**3 / (4 * np.pi**2 * 1.e45)

def deathline_I(P, rho6=1.):
    ''' feed an array of periods and a single value for rho6.
        returns the death line for case I' from Zhang et al. 2000 '''
    return 10**((9./4.)*np.log10(P) - 16.58 + np.log10(rho6))

def deathline_III(P, rho6=1.):
    ''' feed an array of periods and a single value for rho6.
        returns the death line for case III' from Zhang et al. 2000 '''
    return 10**(2*np.log10(P) - 16.52 + np.log10(rho6))

pulsars = np.loadtxt("atnf_catalogue_v1.64_known_Pdot_pulsars.csv", delimiter=",")
pulsars = pulsars.T


magnetars_l = np.loadtxt("McGill_magnetar_catalogue_P_Pdot_limits_noradio.csv", delimiter=",")
magnetars_l = magnetars_l.T
# modify for screen test
magnetars = np.loadtxt("McGill_magnetar_catalogue_P_Pdot_noradio.csv", delimiter=",")
magnetars = magnetars.T

magnetars_r = np.loadtxt("McGill_magnetar_catalogue_P_Pdot_radio.csv", delimiter=",")
magnetars_r = magnetars_r.T

#https://www.nature.com/nature/for-authors/final-submission#:~:text=For%20guidance%2C%20Nature's%20standard%20figure,(120%E2%80%93136%20mm).
fig = plt.figure(figsize=(8.9*cm,8.9*cm))

Pmin = 1.e-3
Pmax = 1.e5

Prange = np.array([Pmin, Pmax])

ax = fig.add_subplot(111)

for B in [1.e10, 1.e12, 1.e14, 1.e16]:
    # de/increase the power of B to change the slope
    # de/increase final divisor to shift left to right
    Pplot = B**0.4 / (10**6.6)
    ax.plot(Prange, Pdot_from_P_B(Prange, B), color="grey", alpha=1.0, linestyle="dotted", zorder=1, lw=0.5)
    ax.text(Pplot, Pdot_from_P_B(Pplot, B), sci_notation(B)+"G", color="black", ha="center", va="top", linespacing=1,  rotation=-30., fontsize=5)#, transform=ax.transAxes)

#labelLines(plt.gca().get_lines(),zorder=2.5)
for E in [1.e26, 1.e30, 1.e34, 1.e38]:
    Pplot = E**(-1./4.) * 10**10
    ax.plot(Prange, Pdot_from_P_E(Prange, E), color="darkcyan", alpha=1.0, linestyle="dotted", zorder=1, lw=0.5)
    ax.text(Pplot, Pdot_from_P_E(Pplot, E), sci_notation(E)+"erg s$^{-1}$", color="darkcyan", ha="center", va="top", linespacing=1,  rotation=60., fontsize=5)#, transform=ax.transAxes)
# FRB death line (Wadiasingh+2020)
#ax.plot([10**-0.7, 10**-0.7], [1.e-10, 1.e-7], color="green", alpha=0.5)
#ax.plot([10**-0.7, 10**1.8], [1.e-10, 1.e-18], color="green", alpha=0.5)

# Pulsar death lines (Zhang et al. 2010)
ax.plot(Prange, deathline_I(Prange, 1.), color="green", alpha=1.0, linestyle="dashed", zorder=0, lw=0.5)
ax.plot(Prange, deathline_III(Prange, 1.), color="green", alpha=1.0, linestyle="dashdot", zorder=0, lw=0.5)
# and another for fun!
ax.scatter(magnetars[0], magnetars[1], facecolors="red", edgecolors="red", marker=".", label="X-ray detected magnetars", s=5, zorder=2)
ax.scatter(magnetars_r[0], magnetars_r[1], facecolors="black", edgecolors="red", marker="o", label="Radio/X-ray magnetars", s=5, zorder=2)
ax.scatter(pulsars[0], pulsars[1], color="black", marker=".", label="Radio pulsars", s=1, zorder=2)
#ax.errorbar(magnetars_l[0], magnetars_l[1], xerr=0, yerr=magnetars_l[1]/2, uplims=1, color="red")
ax.scatter(magnetars_l[0], magnetars_l[1], marker=align_marker(r"$\downarrow$", valign="top"), color="red", s=50, zorder=2)
# Our source
#ax.errorbar(1091.17, 6.e-10, color="blue", xerr=0, yerr=3.e-10, uplims=1, label="GLEAM-X\,J\,162759.5$-$523504.3")
#ax.scatter(1091.17, 6.e-10, marker=align_marker(r"$\downarrow$", valign="top"), color="blue", label="GLEAM-X\,J\,162759.5$-$523504.3", s=150)
# MTP 13
ax.scatter(75.88, 2.25e-13, marker=".", color="black", s=1)
ax.scatter(1091.17, 1.2e-9, marker=align_marker(r"$\downarrow$", valign="top"), color="blue", label="GLEAM-X J162759.5$-$523504.3", s=50)
ax.scatter(1318.57, 8e-14, marker=align_marker(r"$\downarrow$", valign="top"), color="magenta", label="GPM J1839$-$10", s=50)

# Long-period magnetar at the center of RCW 103
#ax.errorbar(6.67*3600, 7.e-10, label="1E 161348–5055", color="magenta", xerr=0, yerr=3.e-10, uplims=1)
ax.scatter(6.67*3600, 7.e-10, marker=align_marker(r"$\downarrow$", valign="top"), color="red", s=50) #, label="1E 161348–5055")
ax.set_xlabel("$P$ / s")
ax.set_ylabel("$\dot{P}$ / s s$^{-1}$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim([1.e-21, 1.e-07])
ax.set_xlim([Pmin, Pmax])
#ax.legend(prop={'size': 8})
ax.legend(frameon=False, handletextpad=0.1, borderaxespad=1, borderpad=0.2, loc='lower right', fontsize=6)

# From Beniamini 2020, I think?
#Pdotgrid = np.loadtxt("Pdotgrid")
#Pgrid = np.loadtxt("Pgrid")
#numdensity = np.loadtxt("numdensity")

# Thanks https://www.tutorialspoint.com/matplotlib/matplotlib_contour_plot.htm
#cvls = ax.contourf(Pgrid, Pdotgrid, np.log10(numdensity), levels = [-3, -2.5, -2, -1.5, -1, -0.5, 0], cmap="bone_r", alpha=0.4, zorder=-10)#, origin="lower", cmap="bone_r", extent=(Pgrid[0],Pdotgrid[0],Pgrid[-1],Pdotgrid[-1]))

# Plot the next two lines if you need to see the levels
#cax = fig.add_subplot(161)
#fig.colorbar(cvls, ax=cax, shrink=0.9)

fig.savefig("P_Pdot.pdf", bbox_inches="tight", dpi=300)
fig.savefig("P_Pdot.png", bbox_inches="tight", dpi=300)
fig.savefig("P_Pdot.eps", bbox_inches="tight", dpi=300)
