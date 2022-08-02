#/usr/bin/env python

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord
import numpy as np
from scipy.signal import convolve2d
from glob import glob
import os
import sys
cm = 1/2.54

try:
    from gleam_x.bin.beam_value_at_radec import beam_value, parse_metafits
    beam_corr = True
except ImportError:
    beam_corr = False

try:
    prefix = sys.argv[1]
except:
    prefix = ""

prepend = "dyn"
debugplot = True

def GetBeamAtCoords(obsid, freq, ra_deg, dec_deg):
#    url = "http://ws.mwatelescope.org/metadata/fits?obs_id=" + str(obs_id)
    t, delays, centfreq, gridnum = parse_metafits(f"{obsid}.metafits")
    beam_x, beam_y = beam_value(ra_deg, dec_deg, t, delays, freq, gridnum,)
    vals = (beam_x + beam_y) / 2
    return vals

def sc(data):
    std = np.std(data)
    std = np.std(data[np.abs(data) < 3*std])
    return std

# https://stackoverflow.com/questions/21921178/binning-a-numpy-array/42024730#42024730
def binArray(data, axis, binstep, binsize, func=np.nanmean):
    data = np.array(data)
    dims = np.array(data.shape)
    argdims = np.arange(data.ndim)
    argdims[0], argdims[axis]= argdims[axis], argdims[0]
    data = data.transpose(argdims)
    data = [func(np.take(data,np.arange(int(i*binstep),int(i*binstep+binsize)),0),0) for i in np.arange(dims[axis]//binstep)]
    data = np.array(data).transpose(argdims)
    return data

hdus = sorted(glob(f"{prefix}_{prepend}-t????-????-image.fits"))
h = hdus[-1]
st = h.split("-t")[1]
tmax = int(st.split("-")[0])
cmax = int(st.split("-")[1])
#    cmax = int(str(hdus[-1])[21:25])
pol=""
#except IndexError:
#    pol="-I"
#    hdus = sorted(glob(f"{prefix}_{prepend}-t????-????{pol}-image.fits"))
#    tmax = int(str(hdus[-1])[16:20])
#    cmax = int(str(hdus[-1])[21:25])

print(tmax, cmax)
ts = np.arange(0, tmax+1)
tsf = [f"{t:04d}" for t in ts]

fs = np.arange(0, cmax+1)
fsf = [f"{f:04d}" for f in fs]

min_nu = fits.open(f"{prefix}_{prepend}-t0000-0000{pol}-image.fits")[0].header["CRVAL3"]-(30.72e6/len(fs)/2)
max_nu = fits.open(f"{prefix}_{prepend}-t0000-{fsf[-1]}{pol}-image.fits")[0].header["CRVAL3"]+(30.72e6/len(fs)/2)

t0 = Time(fits.open(f"{prefix}_{prepend}-t0000-0000{pol}-image.fits")[0].header["DATE-OBS"])
t1 = Time(fits.open(f"{prefix}_{prepend}-t0001-0000{pol}-image.fits")[0].header["DATE-OBS"])
tstep = (t1 - t0).sec

h = fits.open(f"{prefix}_{prepend}-t0000-0000{pol}-image.fits")
x, y = h[0].header["NAXIS1"]/2, h[0].header["NAXIS2"]/2
w = WCS(h[0].header, naxis=2)
coords = w.pixel_to_world(x, y)

nus = []
beam_vals = []

for j in np.arange(0, len(fs)):
    h = fits.open(f"{prefix}_{prepend}-t0000-{fsf[j]}{pol}-image.fits")
    nu = h[0].header["CRVAL3"]
    nus.append(nu)
    beam_vals.append(GetBeamAtCoords(prefix, nu, coords.fk5.ra.value, coords.fk5.dec.value))

np.savetxt(f"{prefix}_frequencies.csv", nus)

vals = []
dates = []

if not os.path.exists(f"{prefix}_{prepend}_dynamic_spectrum.csv"):
# Find maximum
    h = fits.open(hdus[-1])
    box = 5
    arr = np.zeros(shape=(1,1,2*box,2*box), dtype="float32")

    for i in np.arange(0, len(ts)):
        h = fits.open(f"{prefix}_{prepend}-t{tsf[i]}-MFS{pol}-image.fits")
        d = h[0].data
        dates.append(h[0].header["DATE-OBS"])
# Don't include NaN values or it breaks
        if not np.isnan(d).any():
            arr += d[:,:,int(d.shape[2]/2)-box:int(d.shape[2]/2)+box,int(d.shape[3]/2)-box:int(d.shape[3]/2)+box]
# Peak in x, y = location of source
        
#    h[0].data = arr
#    h.writeto("sum.fits")
    peak = np.unravel_index(np.argmax(arr), arr.shape)
# Put it back in the centre
    peak = list(peak)
    if debugplot is True:
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111)
        ax.imshow(np.squeeze(arr).T, aspect='auto', interpolation='none')
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
        ax.scatter(peak[2], peak[3], marker='x', color='red')
        fig.savefig(f"{prefix}_{prepend}_center.png", bbox_inches="tight")
    peak[2] += int(d.shape[2]/2)-box
    peak[3] += int(d.shape[3]/2)-box
    peak = tuple(peak)
    print(peak)
    arr = np.zeros(shape=(len(ts), len(fs)), dtype="float32")
    for i in np.arange(0, len(ts)):
        for j in np.arange(0, len(fs)):
            arr[i, j] = fits.open(f"{prefix}_{prepend}-t{tsf[i]}-{fsf[j]}{pol}-image.fits")[0].data[peak]
        if beam_corr is True:
            arr[i, :] /= beam_vals

    np.savetxt(f"{prefix}_{prepend}_dynamic_spectrum.csv", arr)
    with open(f"{prefix}_{prepend}_timesteps.csv", mode='wt') as myfile:
        myfile.write('\n'.join(dates))
        myfile.write('\n')
else:
    arr = np.loadtxt(f"{prefix}_{prepend}_dynamic_spectrum.csv")

# MeerKAT data cleanup
#m = np.average(np.concatenate([arr[0:100,:],arr[200:300,:]]), axis=0)
#bkg = np.tile(m, (arr.shape[0],1))
#arr -= bkg
        
fig = plt.figure(figsize=(6,3))
ax = fig.add_subplot(111)
#ax.imshow(binArray(arr.T, 1, 2, 2, np.mean))
#kernel_size = (2,2)
#kernel = np.ones(kernel_size)
#ax.imshow(convolve2d(arr.T, kernel, mode='same'))
ax.imshow(arr.T,extent = [0 , tstep*len(ts), max_nu * 1.e-6 , min_nu * 1.e-6 ], aspect='auto', interpolation='none')#, vmin = -0.1, vmax=0.25)
ax.set_xlabel("Time / s")
ax.set_ylabel("Frequency / MHz")
fig.savefig(f"{prefix}_{prepend}_dynspec.pdf", bbox_inches="tight")
