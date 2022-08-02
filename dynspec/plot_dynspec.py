#!/usr/bin/env python

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.time import Time
import numpy as np
from scipy.signal import convolve2d
from glob import glob
import os

import sys

try:
    from gleam_x.bin.beam_value_at_radec import beam_value, parse_metafits
    beam_corr = True
except ImportError:
    beam_corr = False

cm = 1/2.54

try:
    prefix = sys.argv[1]
except:
    prefix = ""

prepend = "pad"


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

try:
    hdus = sorted(glob(f"{prefix}_{prepend}-t????-????-image.fits"))
    tmax = int(str(hdus[-1])[16:20])
    cmax = int(str(hdus[-1])[21:25])
    pol=""
except IndexError:
    pol="-I"
    hdus = sorted(glob(f"{prefix}_{prepend}-t????-????{pol}-image.fits"))
    tmax = int(str(hdus[-1])[16:20])
    cmax = int(str(hdus[-1])[21:25])

print(tmax, cmax)
#ts = np.arange(60*2, 120*2)
ts = np.arange(0, tmax+1)
tsf = [f"{t:04d}" for t in ts]

fs = np.arange(0, cmax+1)
fsf = [f"{f:04d}" for f in fs]

min_nu = fits.open(f"{prefix}_{prepend}-t0000-0000{pol}-image.fits")[0].header["CRVAL3"]-(30.72e6/len(fs)/2)
max_nu = fits.open(f"{prefix}_{prepend}-t0000-{fsf[-1]}{pol}-image.fits")[0].header["CRVAL3"]+(30.72e6/len(fs)/2)

t0 = Time(fits.open(f"{prefix}_{prepend}-t0000-0000{pol}-image.fits")[0].header["DATE-OBS"])
t1 = Time(fits.open(f"{prefix}_{prepend}-t0001-0000{pol}-image.fits")[0].header["DATE-OBS"])
tstep = (t1 - t0).sec

nus = []

for j in np.arange(0, len(fs)):
    h = fits.open(f"{prefix}_{prepend}-t0000-{fsf[j]}{pol}-image.fits")
    nus.append(h[0].header["CRVAL3"])

np.savetxt(f"{prefix}_frequencies.csv", nus)

vals = []
dates = []

if not os.path.exists(f"{prefix}_dynamic_spectrum.csv"):
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
    peak[2] += int(d.shape[2]/2)-box
    peak[3] += int(d.shape[3]/2)-box
    peak = tuple(peak)
    print(peak)
    arr = np.zeros(shape=(len(ts), len(fs)), dtype="float32")
    for i in np.arange(0, len(ts)):
        for j in np.arange(0, len(fs)):
            arr[i, j] = fits.open(f"{prefix}_{prepend}-t{tsf[i]}-{fsf[j]}{pol}-image.fits")[0].data[peak]
            #arr[i, j] = fits.open(f"{prefix}_dyn-t{tsf[i]}-{fsf[j]}{pol}-image.fits")[0].data[:,:,126,125]

    np.savetxt(f"{prefix}_dynamic_spectrum.csv", arr)
    with open(f"{prefix}_timesteps.csv", mode='wt') as myfile:
        myfile.write('\n'.join(dates))
        myfile.write('\n')
else:
    arr = np.loadtxt(f"{prefix}_dynamic_spectrum.csv")
        
fig = plt.figure(figsize=(6,3))
ax = fig.add_subplot(111)
#ax.imshow(binArray(arr.T, 1, 2, 2, np.mean))
#kernel_size = (2,2)
#kernel = np.ones(kernel_size)
#ax.imshow(convolve2d(arr.T, kernel, mode='same'))
ax.imshow(arr.T,extent = [0 , tstep*len(ts), max_nu * 1.e-6 , min_nu * 1.e-6 ], aspect='auto')
ax.set_xlabel("Time / s")
ax.set_ylabel("Frequency / MHz")
fig.savefig(f"{prefix}_dynspec.pdf", bbox_inches="tight")
