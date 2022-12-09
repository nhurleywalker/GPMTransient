#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.time import Time
import matplotlib.cm as cm
from scipy.signal import medfilt
from scipy import optimize
import astropy.constants as const
from astropy.io import ascii

cm = 1/2.54   # centimeters in inches

# Nature requires sans-serif fonts
plt.rcParams.update({
    "text.usetex": False,
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

pol_data = ascii.read("pol_data.csv", format='csv')
#t,ampPeakPIfit,dAmpPeakPIfit,phiPeakPIfit_rm2,dPhiPeakPIfit_rm2,polAngle0Chan_deg,dPolAngle0Chan_deg,polAngleChan_deg,dPolAngleChan_deg
t = Time(pol_data["t"]/(24*60*60), format="mjd")
# to align with observation start time, like other plots in the paper
# have to add 30s of dispersion delay
t0 = Time("2022-07-20T19:11:53", scale="utc", format="isot")

t -= t0
pol_data["t"] = t.value*24*60*60

kwargs = { "lw" : "0.5" , 'markersize':1}
ewargs = { 'markeredgewidth':0.1, 'elinewidth':0.2, 'markersize':1}

# Plot the RMSF
fig = plt.figure(figsize=(8.9*cm, 12.5*cm))
ax1 = fig.add_subplot(411)

plot, = ax1.plot(pol_data["t"], pol_data["ImJy"], marker='.', label="I", color="lightblue", **kwargs)
plot, = ax1.plot(pol_data["t"], 1000.0*pol_data["ampPeakPIfit"], marker='.', color="darkblue", label="P", **kwargs)
ax1.errorbar(pol_data["t"], 1000.0*pol_data["ampPeakPIfit"], yerr=500*pol_data["dAmpPeakPIfit"], fmt="o", color="darkblue", **ewargs)
ax1.legend()
ax1.set_ylabel("$S$ / mJy beam$^{-1}$")

ax2 = fig.add_subplot(412, sharex=ax1)
#plot, = ax2.plot(pol_data["t"], 100000*pol_data["ampPeakPIfit"]/pol_data["ImJy"], marker='.', color="darkblue", label="P%", **kwargs)
# Error budget is dominated by that of the polarised intensity
ratio = 100000*pol_data["ampPeakPIfit"]/pol_data["ImJy"]
yerr = ratio * 0.5*pol_data["dAmpPeakPIfit"]/pol_data["ampPeakPIfit"]
s = np.where(yerr < 50)
ax2.errorbar(pol_data["t"][s], ratio[s], yerr=yerr[s], fmt="o", color="darkblue", **ewargs)
ax2.set_ylabel("$\\frac{P}{I}$ / %")
ax2.set_ylim([0, 100])

ax3 = fig.add_subplot(413, sharex=ax1)
s = np.where(pol_data["dPolAngle0Chan_deg"]<50)
plot, = ax3.plot(pol_data["t"][s], pol_data["polAngle0Chan_deg"][s], marker=".", ls="None", color="darkblue", **kwargs)
plt.errorbar(pol_data["t"][s], pol_data["polAngle0Chan_deg"][s], yerr=pol_data["dPolAngle0Chan_deg"][s], fmt="o", color="darkblue", **ewargs)
ax3.set_ylabel("PA / $^\circ$")

ax4 = fig.add_subplot(414, sharex=ax1)
plot, = ax4.plot(pol_data["t"], pol_data["phiPeakPIfit_rm2"], marker='.', ls="None", color="darkblue",  **kwargs)
plt.errorbar(pol_data["t"], pol_data["phiPeakPIfit_rm2"], yerr=pol_data["dPhiPeakPIfit_rm2"], fmt="o", color="darkblue", **ewargs)
ax4.set_xlabel("Time / s")
ax4.set_ylabel("RM / rad m$^{-2}$")
ax4.set_ylim(-570, -500)

# Since we used a shadarkblue x-axis, we have to be careful when selecting which x-tick labels to hide
# https://stackoverflow.com/questions/4209467/matplotlib-share-x-axis-but-dont-show-x-axis-tick-labels-for-both-just-one#:~:text=This%20is%20a%20common%20gotcha,invisible%20on%20just%20one%20axis.
for ax in [ax1, ax2, ax3]:
    plt.setp(ax.get_xticklabels(), visible=False)
fig.savefig("MeerKAT_corr_pol.pdf", bbox_inches="tight")
plt.close()

