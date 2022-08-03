#!/usr/bin/env python

import numpy as np
from astropy.time import Time
from scipy import interpolate
import matplotlib.pyplot as plt
import csv


mwa200_data = np.loadtxt("1342623496_dyn_dynamic_spectrum.csv")
mwa200_freq = np.loadtxt("1342623496_frequencies.csv")
with open('1342623496_dyn_timesteps.csv', newline='\n') as f:
    reader = csv.reader(f)
    mwa200_time = list(reader)
mwa200_time = np.array([ Time(i[0]).gps for i in mwa200_time ])

mwa88_data = np.loadtxt("1342623792_dyn_dynamic_spectrum.csv").T
mwa88_freq = np.loadtxt("1342623792_frequencies.csv")
with open('1342623792_dyn_timesteps.csv', newline='') as f:
    reader = csv.reader(f)
    mwa88_time = list(reader)
mwa88_time = np.array([ Time(i[0]).gps for i in mwa88_time ])

mwa_delta_freq = mwa88_freq[1] - mwa88_freq[0]
mwa_delta_time = mwa88_time[1] - mwa88_time[0]

print(mwa88_freq[0], mwa200_freq[-1], mwa_delta_freq)
print(mwa200_time[0], mwa88_time[-1], mwa_delta_time)
xbig = np.arange(mwa88_freq[0], mwa200_freq[-1], mwa_delta_freq)
ybig = np.arange(mwa200_time[0], mwa88_time[-1], mwa_delta_time)
#xxbig, yybig = np.meshgrid(xbig, ybig)

arr = np.zeros((len(xbig), len(ybig)))
print(arr.shape)
arr[0:len(mwa200_freq),0:len(mwa200_time)] = np.flipud(mwa200_data.T)
arr[-len(mwa88_freq)-1:-1,-len(mwa88_time)-1:-1] = np.flipud(mwa88_data)

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.imshow(arr, vmin = -0.1, vmax = 10.0, aspect='auto', interpolation='none', extent=[ybig[0]-mwa200_time[0], ybig[-1]-mwa200_time[0], xbig[0]/1.e6, xbig[-1]/1.e6])
ax.set_xlabel("Time / seconds")
ax.set_ylabel("Frequency / MHz")
fig.savefig("split_pulse.pdf", bbox_inches="tight")
