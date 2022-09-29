import sys
import numpy as np
from astropy.io import fits
# A least-squares estimator that gives a useable error estimate
from scipy.optimize import leastsq
from matplotlib import cm as cmap
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

#    "axes.unicode_minus" : False,
#    "pdf.fonttype" : 42, 
#matplotlib.rcParams['ps.fonttype'] = 42
cm = 1/2.54  # centimeters in inches

freq_ref = 200

# Split MWA pulse
dlo = dd.Dynspec(input='../dynspec/1342623792_dyn_dynamic_spectrum.csv', sample_time=0.5, freqlo=72.955, bw=0.08, time_offset=1342623792, transpose=True)
dhi = dd.Dynspec(input='../dynspec/1342623496_dyn_dynamic_spectrum.csv', sample_time=0.5, freqlo=200.955, bw=0.64, time_offset=1342623496, transpose=True)

ylo = dd.parse_yaml('../dynspec/1342623792.yaml')
yhi = dd.parse_yaml('../dynspec/1342623496.yaml')

dlo.freq_ref = freq_ref
dhi.freq_ref = freq_ref
DM = 288 # Not correct, but accounts for scattering at the low end enough that we can sum the profiles for this measurement.

#dlo.dedisperse(DM, freq_ref='low')
dlo.dedisperse(DM)
dhi.dedisperse(DM)

dlo.t = dlo.get_time_at_infinite_frequency()
dhi.t = dhi.get_time_at_infinite_frequency()

fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
fig.suptitle("DM = {:.2f}".format(DM))
dlo.plot(ax[1])
dhi.plot(ax[0])
plt.xlim([dhi.t[0], dlo.t[-1]])
plt.savefig('lineup.png')

# Transpose
# Can only use 72 -- 92 MHz of low-frequency pulse due to time coverage (channels 0 to 250 inclusive)
dat_lo = dlo.dynspec[0:250].T
# The low profile has band gaps which are set to zero, which screws everything up later: set to NaN
dat_lo[dlo.dynspec[0:250].T == 0.0] = np.nan
dat_hi = dhi.dynspec.T

# Likewise for the frequency channels
freq_lo = dlo.f[0:250]
freq_hi = dhi.f

# Helper plot to check that worked (it did)
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.set_xlabel("frequency")
#ax.set_ylabel("flux density (Jy/beam)")
#ax.plot(np.nanmean(dat_lo, axis=0))
#fig.savefig("freqs_lo.png", bbox_inches="tight")
#sys.exit(0)

profile_lo = np.nanmean(dat_lo, axis=1)
# define RMS by a bit of the observation that has no pulse in it
rms_lo = np.nanstd(profile_lo[0:250])
# Define the time range manually so that we integrate over the same amount of the pulse
tstart = 350
tend = 520
mask_1D_lo = np.arange(tstart, tend)
r = np.arange(0,len(profile_lo))

# Helper plot to show what data is therefore included
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("timestep")
ax.set_ylabel("flux density (Jy/beam)")
ax.fill_between(r, -rms_lo, rms_lo, color="red", alpha=0.3)
ax.plot(r, profile_lo)
ax.plot(r[mask_1D_lo],profile_lo[mask_1D_lo], color="green")
fig.savefig("profile_lo.png", bbox_inches="tight")

profile_hi = np.nanmean(dat_hi, axis=1)
rms_hi = np.nanstd(profile_hi[0:250])
# Define the time range manually so that we integrate over the same amount of the pulse
mask_1D_hi = np.arange(tstart, tend)
r = np.arange(0,len(profile_hi))

# Helper plot to show what data is therefore included
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("timestep")
ax.set_ylabel("flux density (Jy/beam)")
ax.fill_between(r, -rms_hi, rms_hi, color="red", alpha=0.3)
ax.plot(r, profile_hi)
ax.plot(r[mask_1D_hi], profile_hi[mask_1D_hi], color="green")
fig.savefig("profile_hi.png", bbox_inches="tight")

# Export the sum of the data over various frequency bands
# For low-frequency data, two 10-MHz chunks
total_lo = np.sum(profile_lo[mask_1D_lo])
total_hi = np.sum(profile_hi[mask_1D_hi])
fluxes = []
fluxes.append(np.nansum(dat_lo[tstart:tend,0:125])/((tend - tstart) * 125))
fluxes.append(np.nansum(dat_lo[tstart:tend,125:250])/((tend - tstart) * 125))

# For the high-frequency data, four 7.68-MHz chunks
fluxes.append(np.nansum(dat_hi[tstart:tend,0:12])/((tend - tstart) * 12))
fluxes.append(np.nansum(dat_hi[tstart:tend,12:24])/((tend - tstart) * 12))
fluxes.append(np.nansum(dat_hi[tstart:tend,24:36])/((tend - tstart) * 12))
fluxes.append(np.nansum(dat_hi[tstart:tend,36:48])/((tend - tstart) * 12))
fluxes = np.array(fluxes)

freqs = np.array([dlo.f[62], dlo.f[187], dhi.f[6], dhi.f[18], dhi.f[30], dhi.f[42]])

# For now, 10% errors
fluxerrs = 0.1*fluxes

ind = np.arange(0,6)

np.savetxt("radio_SED.csv", np.vstack([ind, freqs, fluxes, fluxerrs]).T, delimiter=",", header="ind,freq,flux,fluxerr", fmt=("%0.0d","%3.3f","%5.3f","%5.3f"))

