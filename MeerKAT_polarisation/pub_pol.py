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

#pol_data = ascii.read("pol_data_with_alpha.csv", format='csv')
pol_data = ascii.read("pol_data_after_dedispersion.csv", format='csv')
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
ax1 = fig.add_subplot(511)

plot, = ax1.plot(pol_data["t"], 1000.0*pol_data["IJy"], marker='.', label="I", color="lightblue", **kwargs)
plot, = ax1.plot(pol_data["t"], 1000.0*pol_data["ampPeakPIfit"], marker='.', color="darkblue", label="P", **kwargs)
ax1.errorbar(pol_data["t"], 1000.0*pol_data["ampPeakPIfit"], yerr=1000*pol_data["dAmpPeakPIfit"], fmt="o", color="darkblue", **ewargs)
ax1.legend()
ax1.set_ylabel("$S$ / mJy beam$^{-1}$")

ax2 = fig.add_subplot(512, sharex=ax1)
#plot, = ax2.plot(pol_data["t"], 100000*pol_data["ampPeakPIfit"]/pol_data["ImJy"], marker='.', color="darkblue", label="P%", **kwargs)
# Error budget is dominated by that of the polarised intensity
ratio = 100*pol_data["ampPeakPIfit"]/pol_data["IJy"]
yerr = ratio * pol_data["dAmpPeakPIfit"]/pol_data["ampPeakPIfit"]
s = np.where(yerr < 50)
ax2.errorbar(pol_data["t"][s], ratio[s], yerr=yerr[s], fmt="o", color="darkblue", **ewargs)
ax2.set_ylabel("$\\frac{P}{I}$ / %")
ax2.set_ylim([0, 100])

ax3 = fig.add_subplot(513, sharex=ax1)
#s = np.where(pol_data["dPolAngle0Chan_deg"]<50)
s = np.where(ratio > 0.)
#s = np.where(1000.0*pol_data["ampPeakPIfit"] > 15)
d = np.intersect1d(np.where(ratio < 15.), np.where(1000.0*pol_data["ampPeakPIfit"] > 15))
plot, = ax3.plot(pol_data["t"][s], pol_data["polAngle0Chan_deg"][s], marker=".", ls="None", color="darkblue", **kwargs)
ax3.errorbar(pol_data["t"][s], pol_data["polAngle0Chan_deg"][s], yerr=pol_data["dPolAngle0Chan_deg"][s], fmt="o", color="darkblue", **ewargs)
ax3.scatter(pol_data["t"][d], pol_data["polAngle0Chan_deg"][d], marker="x", color="red", s=10)
ax3.set_ylabel("PA / $^\circ$")

ax4 = fig.add_subplot(514, sharex=ax1)
plot, = ax4.plot(pol_data["t"][s], pol_data["phiPeakPIfit_rm2"][s], marker='.', ls="None", color="darkblue",  **kwargs)
ax4.errorbar(pol_data["t"][s], pol_data["phiPeakPIfit_rm2"][s], yerr=pol_data["dPhiPeakPIfit_rm2"][s], fmt="o", color="darkblue", **ewargs)
ax4.scatter(pol_data["t"][d], pol_data["phiPeakPIfit_rm2"][d], marker="x", color="red", s=10)
ax4.set_xlabel("Time / s")
ax4.set_ylabel("RM / rad m$^{-2}$")
ax4.set_ylim(-570, -500)

ax5 = fig.add_subplot(515, sharex=ax1)
plot, = ax5.plot(pol_data["t"], pol_data["alpha"], marker='.', ls="None", color="darkblue",  **kwargs)
ax5.errorbar(pol_data["t"], pol_data["alpha"], yerr=pol_data["alphaerr"], fmt="o", color="darkblue", **ewargs)
ax5.set_xlabel("Time / s")
ax5.set_ylabel("$\\alpha$")

ax1.set_xlim(100,350)

# Since we used a shared x-axis, we have to be careful when selecting which x-tick labels to hide
# https://stackoverflow.com/questions/4209467/matplotlib-share-x-axis-but-dont-show-x-axis-tick-labels-for-both-just-one#:~:text=This%20is%20a%20common%20gotcha,invisible%20on%20just%20one%20axis.
for ax in [ax1, ax2, ax3, ax4]:
    plt.setp(ax.get_xticklabels(), visible=False)
fig.savefig("MeerKAT_corr_pol.pdf", bbox_inches="tight")
fig.savefig("MeerKAT_corr_pol.png", bbox_inches="tight", dpi=300)
plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(pol_data["ampPeakPIfit"])
ax.set_xlabel("polarised intensity")
ax.set_ylabel("N")
fig.savefig("PI_hist.png", bbox_inches="tight")
