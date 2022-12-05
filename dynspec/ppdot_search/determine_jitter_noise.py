#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt

# Read in the barycentred times of arrival
# Use only the 2022 data
# Work out which pulse is which
# Measure the timing difference between when the pulse arrived and when it "should" have arrived
# Take the root mean square

P = 1318.19578 

# Read in the barycentred times of arrival
arr = np.loadtxt("all_toas_barycentered.tim", skiprows=1, usecols=(1,2))
mjds = arr.T[1]

# Use only the 2022 data
mjds = mjds[mjds > 59770]
# rebaseline
mjds = mjds - mjds[0]

# convert to seconds
toas = mjds * 24 * 60 * 60

# Work out which pulse is which
pulse_nums = np.round(toas / P, decimals=0)

# measure residuals
resid = (toas/P) - pulse_nums

# recentre
resid = resid - np.nanmean(resid)

# Put them in seconds again
resid = P * resid

# measure RMS of residuals
rms = np.std(resid)

# make a helper plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(resid)
ax.axvline(-rms, color='k')
ax.axvline(rms, color='k')
ax.set_xlabel("residual timing errors / s")
fig.savefig("jitter_noise.png")

print("Jitter noise estimate: {rms:2.1f} s")

