#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta
from astropy import units as u
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

v_data = np.loadtxt("MKT_1342378934_V_dynamic_spectrum_dedispersed.csv")
v_profile = np.nanmean(v_data, axis=0)

pol_data = ascii.read("pol_data_dedispersed_with_fits.csv", format='csv')
#t,ampPeakPIfit,dAmpPeakPIfit,phiPeakPIfit_rm2,dPhiPeakPIfit_rm2,polAngle0Chan_deg,dPolAngle0Chan_deg,polAngleChan_deg,dPolAngleChan_deg
t = Time(pol_data["t"]/(24*60*60), format="mjd")
# to align with observation start time, like other plots in the paper
# have to add 30s of dispersion delay
t0 = Time("2022-07-20T19:11:53", scale="utc", format="isot")
t -= t0
pol_data["t"] = t.value*24*60*60

# High-time resolution data doesn't undersample polarisation angle changes
arr = np.loadtxt("t_pa_with_err.txt")
#tp = Time(arr.T[0]/(24*60*60), scale="utc", format="mjd")
tp = TimeDelta(arr.T[0], format='sec')
# According to Yunpeng, the first timestamp is at this MJD
t0 = Time(59780.801915517775, scale="utc", format="mjd")
tp = t0 + tp # Note it is now in MJD, i.e. days not seconds
# Now line up with the rest of the data in the plot
t0 = Time("2022-07-20T19:11:53", scale="utc", format="isot")
tp = tp - t0
# And convert back to seconds
tp = tp.value*24*60*60
# Add an arbitrary scaling to get PA on same scale
offset = 90
pa = arr.T[1] + offset
epa = arr.T[2]



kwargs = { "lw" : "0.5" , 'markersize':1, 'marker' : '.'}
ewargs = { 'markeredgewidth':0.1, 'elinewidth':0.2, 'markersize':1, 'fmt' : "o"}
targs = { "color" : "black" , "fontweight" : "bold"}

# Renormalise the Stokes I fit to zero off-pulse
pol_data["Ifreq0Jy"] -= np.nanmean(pol_data["Ifreq0Jy"][0:50])
rms = np.nanstd(pol_data["Ifreq0Jy"][0:50])

# Plot the RMSF
fig = plt.figure(figsize=(8.9*cm, 12.5*cm))
ax1 = fig.add_subplot(411)
#plot, = ax1.plot(pol_data["t"], 1000.0*pol_data["ImeanJy"], label="ImeanJy", color="blue", **kwargs)
plot, = ax1.plot(pol_data["t"], 1000.0*pol_data["Ifreq0Jy"], label="Stokes I", color="black", **kwargs)
plot, = ax1.plot(pol_data["t"], 1000.0*v_profile[294:568], label="Stokes V", color="magenta", **kwargs)
plot, = ax1.plot(pol_data["t"], 1000.0*pol_data["ampPeakPIfit"], color="blue", label="L", **kwargs)
ax1.errorbar(pol_data["t"], 1000.0*pol_data["ampPeakPIfit"], yerr=1000*pol_data["dAmpPeakPIfit"], color="blue", **ewargs)
ax1.legend(fontsize=5)
ax1.set_ylabel("$S$ (mJy beam$^{-1}$)")

ax2 = fig.add_subplot(412, sharex=ax1)
#plot, = ax2.plot(pol_data["t"], 100000*pol_data["ampPeakPIfit"]/pol_data["ImJy"], color="blue", label="P%", **kwargs)
ratio = 100*pol_data["ampPeakPIfit"]/pol_data["Ifreq0Jy"]
ratio_v = 100*v_profile[294:568]/pol_data["Ifreq0Jy"]
# Errors are quadrature sum of errors on PI and I
yerr = np.abs(ratio * np.sqrt((pol_data["dAmpPeakPIfit"]/pol_data["ampPeakPIfit"])**2 + (rms / pol_data["Ifreq0Jy"])**2))
# Errors on V are quadrature sum of errors on V and I -- which are the same, so sqrt(2) * I
yerr_v = 100*np.abs(np.sqrt(2)*rms / v_profile[294:568])
# Select good polarised intensity points
s = np.where(yerr < 50)
plot, = ax2.plot(pol_data["t"][s], ratio[s], color="blue", ls="None", label="L/I",  **kwargs)
ax2.errorbar(pol_data["t"][s], ratio[s], yerr=yerr[s], color="blue", **ewargs)
# Select good Stokes V intensity points
s = np.where(yerr_v < 50)
plot, = ax2.plot(pol_data["t"][s], ratio_v[s], color="magenta", ls="None", label="V/I", **kwargs)
ax2.errorbar(pol_data["t"][s], ratio_v[s], yerr=yerr_v[s], color="magenta", **ewargs)
ax2.axhline(0, color='k', lw=0.5, alpha=1)
#ax2.errorbar(pol_data["t"], 100*pol_data["fracPol"], yerr=yerr, color="green", **ewargs)
ax2.set_ylabel("frac. pol. (%)")
ax2.set_ylim([-10, 100])
ax2.legend(fontsize=5)

ax3 = fig.add_subplot(413, sharex=ax1)
#s = np.where(pol_data["dPolAngle0Chan_deg"]<50)
#s = np.where(ratio > 0.)
#s = np.where(1000.0*pol_data["ampPeakPIfit"] > 15)
# Low polarisation data where there's enough S/N that we can tell it's depolarised
#d = np.intersect1d(np.where(ratio < 15.), np.where(1000.0*pol_data["ampPeakPIfit"] > 15))
zeropa = 140
ax3.axhline(zeropa+180., color='k', ls='-', lw=0.5, alpha=0.2)
ax3.axhline(zeropa, color='k', ls='-', lw=0.5, alpha=0.2)
ax3.axhline(zeropa+90, color='r', ls=':', lw=0.5, alpha=0.5)
ax3.axhline(zeropa-90, color='r', ls=':', lw=0.5, alpha=0.5)
#plot, = ax3.plot(pol_data["t"][s], pol_data["polAngle0Chan_deg"][s], ls="None", color="blue", **kwargs)
#ax3.errorbar(pol_data["t"][s], pol_data["polAngle0Chan_deg"][s], yerr=pol_data["dPolAngle0Chan_deg"][s], color="blue", **ewargs)
#ax3.scatter(pol_data["t"][d], pol_data["polAngle0Chan_deg"][d], marker="x", color="red", s=10)

# Repeat the plot in the range 180-360 deg
#ax3.plot(pol_data["t"][s], pol_data["polAngle0Chan_deg"][s] + 180, ls="None", color="blue", **kwargs)
#ax3.errorbar(pol_data["t"][s], pol_data["polAngle0Chan_deg"][s] + 180, yerr=pol_data["dPolAngle0Chan_deg"][s], color="blue", **ewargs)
#ax3.scatter(pol_data["t"][d], pol_data["polAngle0Chan_deg"][d] + 180, marker="x", color="red", s=10)

# High-time resolution position angles
s = np.where(pa != offset) 
ax3.plot(tp[s], pa[s], ls="None", color="blue", **kwargs)
ax3.errorbar(tp[s], pa[s], yerr=epa[s], color="blue", **ewargs)
# Repeat the plot in the range 180-360 deg
ax3.plot(tp[s], pa[s] + 180, ls="None", color="blue", **kwargs)
ax3.errorbar(tp[s], pa[s] + 180, yerr=epa[s], color="blue", **ewargs)

ax3.set_xlabel("Time (s)")
ax3.set_ylabel("PA ($^\circ$)")
ax3.set_yticks([0, 180, 360])
ax3.set_ylim([0, 360])

#ax4 = fig.add_subplot(414, sharex=ax1)
#plot, = ax4.plot(pol_data["t"][s], pol_data["phiPeakPIfit_rm2"][s], ls="None", color="blue",  **kwargs)
#ax4.errorbar(pol_data["t"][s], pol_data["phiPeakPIfit_rm2"][s], yerr=pol_data["dPhiPeakPIfit_rm2"][s], color="blue", **ewargs)
#ax4.scatter(pol_data["t"][d], pol_data["phiPeakPIfit_rm2"][d], marker="x", color="red", s=10)
#ax4.set_ylabel("RM / rad m$^{-2}$")
#ax4.set_ylim(-570, -500)

#ax5 = fig.add_subplot(515, sharex=ax1)
#plot, = ax5.plot(pol_data["t"][s], pol_data["alpha"][s], ls="None", color="blue",  **kwargs)
#ax5.errorbar(pol_data["t"][s], pol_data["alpha"][s], yerr=pol_data["alphaerr"][s], fmt="o", color="blue", **ewargs)
##ax5.scatter(pol_data["t"][d], pol_data["alpha"][d], marker="x", color="red", s=10)
#ax5.set_xlabel("Time / s")
#ax5.set_ylabel("$\\alpha$")

# Label a, b... etc as per Nature
x = 105
y = 450
ax1.text(x, y, "a", **targs)
y = 85
ax2.text(x, y, "b", **targs)
y = 310
ax3.text(x, y, "c", **targs)
#y = -510
#ax4.text(x, y, "d", **targs)


# Shared axes so affects all
ax1.set_xlim(100,350)
# Since we used a shared x-axis, we have to be careful when selecting which x-tick labels to hide
# https://stackoverflow.com/questions/4209467/matplotlib-share-x-axis-but-dont-show-x-axis-tick-labels-for-both-just-one#:~:text=This%20is%20a%20common%20gotcha,invisible%20on%20just%20one%20axis.
for ax in [ax1, ax2]:
    plt.setp(ax.get_xticklabels(), visible=False)
fig.savefig("MeerKAT_corr_pol.pdf", bbox_inches="tight")
fig.savefig("MeerKAT_corr_pol.png", bbox_inches="tight", dpi=300)
plt.close()

#pol_data = ascii.read("pol_data_after_dedispersion.csv", format='csv')
#fig = plt.figure()
#ax = fig.add_subplot(111)
#t = np.intersect1d(np.where(pol_data["tindex"]>=363), np.where(pol_data["tindex"]<=465))
#ax.plot(pol_data["IJy"][t], label="after dedispersion")
#pol_data = ascii.read("pol_data_with_alpha.csv", format='csv')
#t = np.intersect1d(np.where(pol_data["tindex"]>=363), np.where(pol_data["tindex"]<=465))
#ax.plot(pol_data["IJy"][t], label="before dedispersion")
#ax.set_ylabel("Stokes I flux density / Jy")
#ax.set_xlabel("timestep")
#ax.legend()
#fig.savefig("I_comparison.png", bbox_inches="tight")
