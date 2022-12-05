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
import pandas as pd
import scipy.stats

from plot_P_Pdot import align_marker

def P(f0):
   return 1/f0

def Pdot(f0, f1):
   return - f1 / (f0**2)

# https://www.cv.nrao.edu/~sransom/web/Ch6.html
def Edot(P, Pdot, I=1.e45):
   return 4 * np.pi**2 * I * Pdot / P**3

def B(P, Pdot):
   return 3.2e19 * np.sqrt(P*Pdot)

def tau(P, Pdot):
   return P / (2*Pdot)

def s_to_Myr(t):
   return t / (60 * 60 * 24 * 365.25 * 1.e6)

# Nature requires sans-serif fonts
plt.rcParams.update({
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

cm = 1/2.54  # centimeters in inches

#df = pd.read_csv("chi2_grid_orig.csv", delimiter=",")
df = pd.read_csv("chi2_grid.csv", delimiter=",")
#df.columns = ["Name", "RA", "Dec"]
print(df.keys())
xsize = len(np.unique(df['F0']))
ysize = len(np.unique(df['F1']))
arr = np.array(df['chi2']).reshape((ysize, xsize))
F0 = np.array(df['F0']).reshape((ysize, xsize))
F1 = np.array(df['F1']).reshape((ysize, xsize))

# First island: minimum chi^2 in main ellipse:
best_f0 = df['F0'][np.argmin(arr)]
best_f1 = df['F1'][np.argmin(arr)]
best_P = P(best_f0)
best_Pdot = Pdot(best_f0, best_f1)
bestfit = np.min(arr)

# Second island: minimum chi^2 in lower-left ellipse:
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
#CIs = (scipy.stats.norm().cdf(nsigma) - 0.5) * 2
#print(f"Confidence intervals for {nsigma} sigma: {CIs}")
# chi^2 random variable for 1 parameters
#rv = scipy.stats.chi2(1)
#contour_levels_1param = rv.ppf(CIs)
#print(f"Contour levels for {nsigma} sigma and 1 parameter: {contour_levels_1param}")

# Plot the grid/contour results
#fig, ax = plt.subplots(figsize=(16, 9))
# Just plot the values offset from the best-fit values



#https://www.nature.com/nature/for-authors/final-submission#:~:text=For%20guidance%2C%20Nature's%20standard%20figure,(120%E2%80%93136%20mm).
fig = plt.figure(figsize=(8.9*cm,7.8*cm))
ax = fig.add_subplot(111)
twod = ax.contour(
    1.e9*(F0 - best_f0),
    1.e18*F1,
    arr - bestfit,
    levels=contour_levels,
    colors="b",
    linewidths=[0.5],
)
fmt = {}
strs = ['$1\\sigma$', '$2\\sigma$', '$3\\sigma$']
for l, s in zip(twod.levels, strs):
    fmt[l] = s
ax.clabel(twod, twod.levels, inline=True, fmt=fmt, fontsize=5)

xy = twod.collections[2].get_paths()[0].vertices
sig3_f0, sig3_f1 = xy[np.argmin(xy.T[1])]
sig3_P = P(sig3_f0/1.e9 + best_f0)
sig3_Pdot = Pdot(sig3_f0/1.e9 + best_f0, sig3_f1/1.e18)
sig3_B = B(sig3_P, sig3_Pdot)
sig3_Edot = Edot(sig3_P, sig3_Pdot)
sig3_tau = s_to_Myr(tau(sig3_P, sig3_Pdot))

xy = twod.collections[1].get_paths()[0].vertices
sig2_f0, sig2_f1 = xy[np.argmin(xy.T[1])]
sig2_P = P(sig2_f0/1.e9 + best_f0)
sig2_Pdot = Pdot(sig2_f0/1.e9 + best_f0, sig2_f1/1.e18)
sig2_B = B(sig2_P, sig2_Pdot)
sig2_Edot = Edot(sig2_P, sig2_Pdot)
sig2_tau = s_to_Myr(tau(sig2_P, sig2_Pdot))

xy = twod.collections[0].get_paths()[0].vertices
sig1_f0, sig1_f1 = xy[np.argmin(xy.T[1])]
sig1_P = P(sig1_f0/1.e9 + best_f0)
sig1_Pdot = Pdot(sig1_f0/1.e9 + best_f0, sig1_f1/1.e18)
sig1_B = B(sig1_P, sig1_Pdot)
sig1_Edot = Edot(sig1_P, sig1_Pdot)
sig1_tau = s_to_Myr(tau(sig1_P, sig1_Pdot))

#print(np.nanmin(df['F0'])- best_f0)
#print(np.nanmax(df['F0'])- best_f0)
im = ax.imshow(np.log10(arr), origin="lower", extent=[1.e9*(np.nanmin(df['F0'])- best_f0), 1.e9*(np.nanmax(df['F0'])-best_f0), 1.e18*np.nanmin(df['F1']), 1.e18*np.nanmax(df['F1'])], interpolation="none", aspect="auto", cmap="plasma_r")
ax.set_xlabel(r'$\Delta f$  / nHz')
ax.set_ylabel(r'$\Delta \dot{f}$ / $10^{-18}$')
ax.scatter(best_f0 - best_f0, 1.e18*best_f1, marker="+", zorder=30, lw=0.5, s=5)
ax.scatter(sig1_f0 - best_f0/1.e9, sig1_f1, marker=align_marker(r"$\uparrow$", valign="bottom"), color="magenta", alpha=1, zorder=30, lw=0.45)
ax.scatter(sig2_f0 - best_f0/1.e9, sig2_f1, marker=align_marker(r"$\uparrow$", valign="bottom"), color="magenta", alpha=0.5, zorder=30, lw=0.45)
ax.scatter(sig3_f0 - best_f0/1.e9, sig3_f1, marker=align_marker(r"$\uparrow$", valign="bottom"), color="magenta", alpha=0.2, zorder=30, lw=0.45)
plt.colorbar(im, label=r"$\log_{10}(\chi^2)$")
fig.savefig("Ppdot_search.pdf", bbox_inches="tight", dpi=300)

print(f"Best F0 = {best_f0}, Best F1 = {best_f1}")
print(f"Best P = {best_P}, Best Pdot = {best_Pdot}")

print(f"1-sigma limit F0 = {sig1_f0/1.e9 + best_f0}, 1-sigma limit F1 = {sig1_f1/1.e18}")
print(f"1-sigma limit P = {sig1_P}, 1-sigma limit Pdot = {sig1_Pdot}")
print(f"1-sigma limit Edot = {sig1_Edot:2.2g} erg/s, B = {sig1_B:2.2g} G, tau = {sig1_tau:2.2g} Myr")

print(f"2-sigma limit F0 = {sig2_f0/1.e9 + best_f0}, 2-sigma limit F1 = {sig2_f1/1.e18}")
print(f"2-sigma limit P = {sig2_P}, 2-sigma limit Pdot = {sig2_Pdot}")
print(f"2-sigma limit Edot = {sig2_Edot:2.2g} erg/s, B = {sig2_B:2.2g} G, tau = {sig2_tau:2.2g} Myr")

print(f"3-sigma limit F0 = {sig3_f0/1.e9 + best_f0}, 3-sigma limit F1 = {sig3_f1/1.e18}")
print(f"3-sigma limit P = {sig3_P}, 3-sigma limit Pdot = {sig3_Pdot}")
print(f"3-sigma limit Edot = {sig3_Edot:2.2g} erg/s, B = {sig3_B:2.2g} G, tau = {sig3_tau:2.2g} Myr")



#fig.savefig("P_Pdot.pdf", bbox_inches="tight", dpi=300)

#fig.savefig("P_Pdot.pdf", bbox_inches="tight", dpi=300)

#fig.savefig("P_Pdot.pdf", bbox_inches="tight", dpi=300)
#fig.savefig("P_Pdot.png", bbox_inches="tight", dpi=300)
#fig.savefig("P_Pdot.eps", bbox_inches="tight", dpi=300)
