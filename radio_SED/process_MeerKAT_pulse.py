import sys
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
# Make sure to do export PYTHONPATH=$HOME/Programs/GPMTransient/dynspec/
import dedisperse_dynspec as dd

# Nature requires sans-serif fonts
plt.rcParams.update({
    "text.usetex": False,
    "font.size": 7,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

cm = 1/2.54  # centimeters in inches

arr = np.loadtxt("../dynspec/MKT_1342379534_dyn_dynamic_spectrum.csv")
# MeerKAT data cleanup
m = np.average(np.concatenate([arr[0:100,:],arr[200:300,:]]), axis=0)
bkg = np.tile(m, (arr.shape[0],1))
arr -= bkg

d = dd.Dynspec(dynspec=arr, sample_time=2, freqlo=568.171875, bw=6.375, time_offset=1342379534, transpose=True)

d.freq_ref = 1000
DM = 275

d.dedisperse(DM)

d.t = d.get_time_at_infinite_frequency()

fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True)
fig.suptitle("DM = {:.2f}".format(DM))
# Some manual flagging
d.dynspec[0:3] = np.nan
d.dynspec[56:60] = np.nan
d.plot(ax)
# So much great S/N... we have 76 channels, although the bottom 3 and top one are flagged -> 72 channels
# -> 8 x 9channel-chunks
for i in range(0, 9):
    ax.axhline(d.f[3+(9*i)], color='k', lw=0.5)
plt.savefig('MeerKAT.png')

# Transpose
dat = d.dynspec.T

profile = np.nanmean(dat, axis=1)
# define RMS by a bit of the observation that has no pulse in it
rms = np.nanstd(profile[200:300])
# Define the time range manually so that we integrate over the same amount of the pulse
tstart = 53
tend = 150
mask_1D = np.arange(tstart, tend)
r = np.arange(0,len(profile))

# Helper plot to show what data is therefore included
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("timestep")
ax.set_ylabel("flux density (Jy/beam)")
ax.fill_between(r, -rms, rms, color="red", alpha=0.3)
ax.plot(r, profile)
ax.plot(r[mask_1D],profile[mask_1D], color="green")
fig.savefig("profile.png", bbox_inches="tight")

fluxes = []
freqs = []
for i in range(0, 8):
    fluxes.append(np.nansum(dat[tstart:tend,3+(9*i):3+(9*(i+1))])/((tend - tstart) * 9))
    freqs.append(d.f[int(3+(9*(i+0.5)))])

freqs = np.array(freqs)
fluxes = np.array(fluxes)

# For now, 10% errors
fluxerrs = 0.1*fluxes

ind = np.arange(0, 8)

np.savetxt("radio_SED3.csv", np.vstack([ind, freqs, fluxes, fluxerrs]).T, delimiter=",", header="ind,freq,flux,fluxerr", fmt=("%0.0d","%3.3f","%5.3f","%5.3f"))

