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

def P(f0):
   return 1/f0

def Pdot(f0, f1):
   return - f1 / (f0**2)

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
fig = plt.figure(figsize=(8.9*cm,8.9*cm))
ax = fig.add_subplot(111)
twod = ax.contour(
    F0,
    F1,
    arr - bestfit,
    levels=contour_levels,
    colors="b",
    linewidths=[0.5],
)

xy = twod.collections[2].get_paths()[0].vertices
sig3_f0, sig3_f1 = xy[np.argmin(xy.T[1])]
sig3_P = P(sig3_f0)
sig3_Pdot = Pdot(sig3_f0, sig3_f1)
xy = twod.collections[1].get_paths()[0].vertices
sig2_f0, sig2_f1 = xy[np.argmin(xy.T[1])]
sig2_P = P(sig2_f0)
sig2_Pdot = Pdot(sig2_f0, sig2_f1)
xy = twod.collections[0].get_paths()[0].vertices
sig1_f0, sig1_f1 = xy[np.argmin(xy.T[1])]
sig1_P = P(sig1_f0)
sig1_Pdot = Pdot(sig1_f0, sig1_f1)

im = ax.imshow(np.log(arr), origin="lower", extent=[np.nanmin(df['F0']), np.nanmax(df['F0']), np.nanmin(df['F1']), np.nanmax(df['F1'])], interpolation="none", aspect="auto", cmap="plasma_r")
#ax.set_xlim([2.2e-9 + 7.5861e-4, 3.3e-9 + 7.5861e-4])
#ax.set_ylim([-0.6e-17, 0.6e-17])
ax.set_xlabel('F0')
ax.set_ylabel('F1')
ax.scatter(best_f0, best_f1, marker="+", zorder=30, lw=0.5)
ax.scatter(sig1_f0, sig1_f1, marker="x", zorder=30, lw=0.5)
ax.scatter(sig2_f0, sig2_f1, marker="x", zorder=30, lw=0.5)
ax.scatter(sig3_f0, sig3_f1, marker="x", zorder=30, lw=0.5)
plt.colorbar(im)
fig.savefig("test.png", bbox_inches="tight", dpi=300)

print(f"Best F0 = {best_f0}, Best F1 = {best_f1}")
print(f"Best P = {best_P}, Best Pdot = {best_Pdot}")

print(f"1-sigma limit F0 = {sig1_f0}, 1-sigma limit F1 = {sig1_f1}")
print(f"1-sigma limit P = {sig1_P}, 1-sigma limit Pdot = {sig1_Pdot}")

print(f"2-sigma limit F0 = {sig2_f0}, 2-sigma limit F1 = {sig2_f1}")
print(f"2-sigma limit P = {sig2_P}, 2-sigma limit Pdot = {sig2_Pdot}")

print(f"3-sigma limit F0 = {sig3_f0}, 3-sigma limit F1 = {sig3_f1}")
print(f"3-sigma limit P = {sig3_P}, 3-sigma limit Pdot = {sig3_Pdot}")




#fig.savefig("P_Pdot.pdf", bbox_inches="tight", dpi=300)

#fig.savefig("P_Pdot.pdf", bbox_inches="tight", dpi=300)

#fig.savefig("P_Pdot.pdf", bbox_inches="tight", dpi=300)
#fig.savefig("P_Pdot.png", bbox_inches="tight", dpi=300)
#fig.savefig("P_Pdot.eps", bbox_inches="tight", dpi=300)
