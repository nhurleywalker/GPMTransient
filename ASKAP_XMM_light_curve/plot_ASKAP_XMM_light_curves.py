#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
import pandas as pd
import math

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as font_manager
import matplotlib.gridspec as gridspec

from numpy import arange

from scipy import interpolate
from scipy import optimize
from scipy.integrate import quad
import scipy.special

import os
from matplotlib import patches
from matplotlib import colors
from matplotlib.lines import Line2D
from matplotlib import rc, rcParams
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

# Nature requires sans-serif fonts
plt.rcParams.update({
    "text.usetex": False,
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

single_col = 8.9 # cm
double_col = 17.8 # cm
def cm2inch(value):
    return value/2.54

# For convenient plotting without long MJDs
offset = 59836

# Read in the data
askap = np.genfromtxt('ASKAP_lightcurve_barycentred.txt')
askap_time = askap[:,0] - offset
askap_flux = askap[:,1]

xmm = np.genfromtxt('gpm1839_xmm-EPIC_lcurve_30min.txt')
xmm_time = xmm[:,0] - offset
xmm_time_err = xmm[:,1]
xmm_flux = xmm[:,2]
xmm_flux_err = xmm[:,3]

# Using gridspec to define the figure, to remove the vertical white space between the 3 panels. 
#fig = plt.figure(figsize=(20,9))
fig = plt.figure(figsize=(cm2inch(single_col),cm2inch(single_col)))
gs1 = gridspec.GridSpec(2, 1, height_ratios=[1,1])
gs1.update(wspace=0.0, hspace=0.0) # removing the white space

ax = plt.subplot(gs1[0])

ax.plot(askap_time, askap_flux, lw=0.5, color = 'darkblue', label = 'Radio flux density (mJy)')

ax.xaxis.set_ticks_position('both')
ax.set_xticklabels([])
ax.get_yaxis().set_tick_params(direction='in', which='both')
ax.get_xaxis().set_tick_params(direction='in', which='both', top=False)

ax.set_yscale('linear')
#ax.set_xlim(np.min(askap_time), np.max(askap_time))
# To include both the ASKAP and XMM overlapping data
ax.set_xlim(0.4343, 0.4973)
ax.set_ylim(-14.0, 145.0)
#ax.yaxis.set_major_locator(plt.FixedLocator([0, 20, 40, 60, 80, 100, 120, 140]))

ax.set_ylabel('Radio flux density (mJy)')#, labelpad=55)

bx = plt.subplot(gs1[1], sharex=ax)

bx.errorbar(xmm_time ,xmm_flux*1.e3,xerr=xmm_time_err, yerr=xmm_flux_err*1.e3, fmt = 'o', lw=0.5, markersize=1,
            mec = 'darkred', color = 'darkred', ecolor='darkred', capsize=3, label = 'XMM-Newton/EPIC (0.3-10 keV)')

plt.axhline(y=0.002*1.e3, color='darkblue', linestyle='dotted', lw=0.5)

#bx.tick_params(width=3, length=12, axis='both', which='major', pad=5)
#bx.tick_params(length=7, width=2, axis='both', which='minor', pad=5)

bx.xaxis.set_ticks_position('both')
bx.get_yaxis().set_tick_params(direction='in', which='both')
bx.get_xaxis().set_tick_params(direction='in', which='both', top=False)

bx.set_yscale('linear')
bx.set_ylim(-0.005*1.e3, 0.018*1.e3)
#bx.yaxis.set_major_locator(plt.FixedLocator([0.0, 0.005, 0.010, 0.015]))

bx.set_ylabel('XMM-Newton/EPIC \n net count rate (10$^{-3}$ counts s$^{-1}$)')#, labelpad=10)
bx.set_xlabel(r'Epoch (BMJD$_{{\rm TDB}}$ - {0})'.format(offset)) # Bottom panel: defining the x-axis label

#bx.ticklabel_format(useOffset=False)
#bx.xaxis.set_major_locator(plt.FixedLocator([59836.44, 59836.45, 59836.46, 59836.47, 59836.48, 59836.49]))

##################################

    
####################################################################################################################
# Saving the figure.
plt.tight_layout()
plt.savefig('GPM1839_askap_xmm_lcurves.pdf')


# In[ ]:




