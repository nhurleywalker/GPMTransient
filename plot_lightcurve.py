#!/usr/bin/env python

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.time import Time
import numpy as np
from glob import glob

import sys

try:
    prefix = sys.argv[1]
except:
    prefix = ""

def sc(data):
    std = np.std(data)
    std = np.std(data[np.abs(data) < 3*std])
    return std


fig = plt.figure()
ax = fig.add_subplot(111)
#for pol, alpha in zip(["I", "Q", "U", "V"], [1, 0.5, 0.5, 0.5]):
pol="I"
alpha=1
#hdus = sorted(glob("*-t05??-image.fits"))
hdus = sorted(glob(f"{prefix}_time-t????-image.fits"))

val = []
std = []
#epo = []
#freqs = []

# Get time delta in seconds
d0 = Time(fits.open(hdus[0])[0].header["DATE-OBS"])
d1 = Time(fits.open(hdus[1])[0].header["DATE-OBS"])
ts = int((d1 - d0).sec)

for hdu in hdus:
    h = fits.open(hdu)
    try:
        val.append(np.average(h[0].data[:,:,124:127,124:127]))
#        epo.append(int(hdu[0:10]))
        std.append(sc(h[0].data))
#        freqs.append(h[0].header["CRVAL3"]/1.e6)
    except ValueError:
        pass

#scat = ax.scatter(epo, val)#, c = freqs)
#clb = plt.colorbar(scat, label="Observing frequency / MHz")
ax.errorbar(ts*np.arange(0,len(val)), val, yerr=std, label=pol, alpha=0.5, lw=0.5, color='k')#, marker='', ls='', zorder=0)
#ax.scatter(ts*np.arange(0,len(val)), val, alpha=alpha, marker='.', zorder=0, color='k')
ax.plot(ts*np.arange(0,len(val)), val, alpha=alpha, lw=0.1, zorder=0, color='k')
ax.set_xlabel("Time (s)")
ax.set_ylabel("Brightness (Jy/beam)")
ax.set_ylim([-0.2,1.0])
#    ax.legend()
fig.savefig(f"{prefix}_light_curve.pdf", bbox_inches="tight")
fig.savefig(f"{prefix}_light_curve.png", bbox_inches="tight")
