#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.time import Time
#import matplotlib.cm as cm
from scipy.signal import medfilt
from scipy import optimize
import astropy.constants as const
import xarray
from astropy.io import ascii

cm = 1/2.54   # centimeters in inches

# Nature requires sans-serif fonts
plt.rcParams.update({
    "text.usetex": False,
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

ix = xarray.open_dataset("1658342033_sdp_l0_GPM1839-10_polcal_scan7.ms_StokesI-kata.nc")

t = Time("2022-07-20T19:11:53", scale="utc", format="isot")

t0 = int(Time(ix["time"][0].data/(24*60*60), format="mjd").gps)
t1 = int(Time(ix["time"][-1].data/(24*60*60), format="mjd").gps)
nu0 = ix["freq"][0].data/1.e6
nu1 = ix["freq"][-1].data/1.e6

I = np.loadtxt("MKT_1342378934_I_dynamic_spectrum_dedispersed.csv")
Q = np.loadtxt("MKT_1342378934_Q_dynamic_spectrum_dedispersed.csv")
U = np.loadtxt("MKT_1342378934_U_dynamic_spectrum_dedispersed.csv")
V = np.loadtxt("MKT_1342378934_V_dynamic_spectrum_dedispersed.csv")

iargs = { "origin" : "lower", "vmin" : -0.1, "vmax" : 0.5, "extent" : [t0-t.gps, t1-t.gps, nu0, nu1], "aspect" : 'auto', "interpolation" : 'none'}
pargs = { "origin" : "lower", "vmin" : -0.1, "vmax" : 0.1, "extent" : [t0-t.gps, t1-t.gps, nu0, nu1], "aspect" : 'auto', "interpolation" : 'none'}
targs = { "color" : "white" , "fontweight" : "bold"}

fig = plt.figure(figsize=(17.8*cm, 17.8*cm))
x , y = 105, 995
ax1 = fig.add_subplot(221)
ax1.imshow(I, **iargs)
ax1.text(x, y, "a", **targs)
ax2 = fig.add_subplot(222)
ax2.imshow(V, **pargs)
ax2.text(x, y, "b", **targs)
ax3 = fig.add_subplot(223)
ax3.imshow(Q, **pargs)
ax3.text(x, y, "c", **targs)
ax4 = fig.add_subplot(224)
ax4.imshow(U, **pargs)
ax4.text(x, y, "d", **targs)

for ax in [ax1, ax2, ax3, ax4]:
# First pulse
    ax.set_xlim(100, 350)
# Second pulse
#    ax.set_xlim(1400, 1500)
    ax.set_ylim(590, 1015)

plt.subplots_adjust(hspace = 0, wspace = 0)

for ax in [ax1, ax2]:
    plt.setp(ax.get_xticklabels(), visible=False)
for ax in [ax3, ax4]:
    ax.set_xlabel("Time / s")
for ax in [ax2, ax4]:
    plt.setp(ax.get_yticklabels(), visible=False)
for ax in [ax1, ax3]:
    ax.set_ylabel("Frequency / MHz")

# Prune overlapping tick labels

ax3.xaxis.set_major_locator(MaxNLocator(prune='upper'))
ax4.xaxis.set_major_locator(MaxNLocator(prune='lower'))
#plt.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))

fig.savefig("MeerKAT_IQUV.pdf", bbox_inches="tight", dpi=300)

