#!/usr/bin/env python

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
import numpy as np
from glob import glob

import sys

try:
    from gleam_x.bin.beam_value_at_radec import beam_value, parse_metafits
    beam_corr = True
except ImportError:
    beam_corr = False

try:
    prefix = sys.argv[1]
except:
    prefix = ""

def sc(data):
    std = np.std(data)
    std = np.std(data[np.abs(data) < 3*std])
    return std

def GetBeamAtCoords(obsid, freq, ra_deg, dec_deg):
#    url = "http://ws.mwatelescope.org/metadata/fits?obs_id=" + str(obs_id)
    t, delays, centfreq, gridnum = parse_metafits(f"{obsid}.metafits")
    beam_x, beam_y = beam_value(ra_deg, dec_deg, t, delays, freq, gridnum,)
    vals = (beam_x + beam_y) / 2
    return vals

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
ts = round((d1 - d0).sec)

# Find maximum
h = fits.open(hdus[0])
box = 1
arr = np.zeros(shape=(1,1,2*box,2*box), dtype="float32")

# Find frequency
nu = h[0].header["CRVAL3"]

# And WCS
w = WCS(h[0].header, naxis=2)

# Find location of peak pixel
for hdu in hdus:
    h = fits.open(hdu)
    d = h[0].data
# Don't include NaN values or it breaks
    if not np.isnan(d).any():
        arr += d[:,:,int(d.shape[2]/2)-box:int(d.shape[2]/2)+box,int(d.shape[3]/2)-box:int(d.shape[3]/2)+box]
# Peak in x, y = location of source
peak = np.unravel_index(np.argmax(arr), arr.shape)
# Put it back in the centre
peak = list(peak)
peak[2] += int(d.shape[2]/2)-box
peak[3] += int(d.shape[3]/2)-box
coords = w.pixel_to_world(peak[3],peak[2])
peak = tuple(peak)

beam_vals = []
for hdu in hdus:
    h = fits.open(hdu)
    try:
        val.append(h[0].data[peak])
#        epo.append(int(hdu[0:10]))
        std.append(sc(h[0].data))
#        freqs.append(h[0].header["CRVAL3"]/1.e6)
        if beam_corr is True:
            beam_vals.append(GetBeamAtCoords(prefix, nu, coords.fk5.ra.value, coords.fk5.dec.value ))
        else:
            beam_vals.append(1)

    except ValueError:
        pass

if beam_corr is True:
    val = np.array(val) / np.array(beam_vals)
    std = np.array(std) / np.array(beam_vals)

#scat = ax.scatter(epo, val)#, c = freqs)
#clb = plt.colorbar(scat, label="Observing frequency / MHz")
ax.errorbar(ts*np.arange(0,len(val)), val, yerr=std, label=pol, alpha=0.5, lw=0.5, color='k')#, marker='', ls='', zorder=0)
#ax.scatter(ts*np.arange(0,len(val)), val, alpha=alpha, marker='.', zorder=0, color='k')
ax.plot(ts*np.arange(0,len(val)), val, alpha=alpha, lw=0.1, zorder=0, color='k')
ax.set_xlabel("Time (s)")
ax.set_ylabel("Brightness (Jy/beam)")
#ax.set_ylim([-0.2,1.0])
#    ax.legend()
fig.savefig(f"{prefix}_light_curve.pdf", bbox_inches="tight")
fig.savefig(f"{prefix}_light_curve.png", bbox_inches="tight")

time_arr = np.arange(d0.gps, d0.gps + len(val)*ts, ts)
np.savetxt(f"{prefix}_light_curve.csv", np.vstack((time_arr, val)).T, fmt="%0.0d %f")
