#!/usr/bin/env python

import xarray
import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.time import Time
import matplotlib.cm as cm
from scipy.signal import medfilt
from scipy import optimize
import astropy.constants as const

def getFDF2D(dataQ, dataU, freqs, startPhi, stopPhi, dPhi, dType='float32'):
    nX = dataQ.shape[0]
    # Calculate the RM sampling
    phiArr = np.arange(startPhi, stopPhi, dPhi)
    nPhi = phiArr.shape[0]

    # Calculate the frequency and lambda sampling
    lamSqArr = np.power(const.c.value / np.array(freqs), 2.0)

    # Calculate the dimensions of the output RM cube
    nPhi = len(phiArr)

    # Initialise the complex Faraday Dispersion Function (FDF)
#    FDF = np.ndarray((nPhi), dtype='complex')
    FDFcube = np.ndarray((nX, nPhi), dtype='complex')

    # Assume uniform weighting
    wtArr = np.ones(len(lamSqArr), dtype=dType)

    K = 1.0 / np.nansum(wtArr)

    # Get the weighted mean of the LambdaSq distribution (B&dB Eqn. 32)
    lam0Sq = K * np.nansum(wtArr * lamSqArr)

    # Mininize the number of inner-loop operations by calculating the
    # argument of the EXP term in B&dB Eqns. (25) and (36) for the FDF
    a = (-2.0 * 1.0j * phiArr)
    b = (lamSqArr - lam0Sq) 
    arg = np.exp( np.outer(a, b) )

    # Create a weighted complex polarised surface-brightness cube
    # i.e., observed polarised surface brightness, B&dB Eqns. (8) and (14)
    PobsCube = (np.array(dataQ) + 1.0j * np.array(dataU)) * wtArr

    for i in range(nX):
        # Calculate the Faraday Dispersion Function
        # B&dB Eqns. (25) and (36)
        FDFcube[i,:] = K * np.nansum(PobsCube[i,:] * arg, 1)
        
    return FDFcube, phiArr, lam0Sq

def get_fdf_2D(times, freqs, It, Qt, Ut, Vt, startPhi, dPhi):
    stopPhi = -startPhi+dPhi

    return getFDF2D(Qt, Ut, freqs, startPhi, stopPhi, dPhi)

# Plot the Faraday Dispersion Function for the given time integration
def plot_fdf_2D(t, freqs, I, Q, U, V, startPhi, dPhi):
    averI = np.nanmean(I, axis=1)  # Average DS in frequency and convert to mJy
    averV = np.nanmean(V, axis=1)
    fdf, phi, lam0Sq = get_fdf_2D(t, freqs, I, Q, U, V, startPhi, dPhi)
    print(fdf.shape)
    sigma = np.std(fdf[:200][0]) # rough estimate of noise
    pi = np.max(np.abs(fdf), axis=1)
    madmax = np.where(pi<8.5 * sigma)
    pi[madmax] = np.nan
    pimaxpos = np.argmax(np.abs(fdf), axis=1)
    phi_peak = phi[pimaxpos]
    polAngleChan_deg = np.zeros(t.shape, dtype=np.float)
    for index in range(len(t)):
#        if pi[index] > 0.005:
        polAngleChan_deg[index] = 0.5 * np.degrees(np.arctan2(fdf[index,pimaxpos[index]].imag, fdf[index,pimaxpos[index]].real)) % 180
        
    polAngle0Chan_deg = (np.degrees(np.radians(polAngleChan_deg) - phi_peak * lam0Sq)) % 180
    # Error propagation for arctan
    # d(arctan(x/y))/dx = y / (x^2 + y^2)
    # So error on arctan(imag / real) is ((real*err_real + imag*err_imag) / (imag^2 + real^2))) ?
    # 
    
    polAngle0Chan_deg[madmax] = np.nan
    phi_peak[madmax] = np.nan
    # Plot the RMSF
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    plot, = ax1.plot(t, polAngle0Chan_deg, marker=".", ls="None", color="grey", label="PA")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("PA (deg)")
    #ax1.set_xlim(5165061225, 5165061425)
    ax1.set_xlim(t[294], t[568])
    ax1.set_ylim(140.0, 180.0)
    
#    ax1.set_xlim(t[0], t[-1])
    ax1 = fig.add_subplot(312)
    plot, = ax1.plot(t, 1000.0*pi, marker='.', color="grey", label="PI")
    plot, = ax1.plot(t, 1000.0*averI, marker='.', color="black", label="I")
    plot, = ax1.plot(t, 1000.0*averV, marker='.', color="green", label="V")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Flux (mJy beam$^{-1}$)")
    #ax1.set_xlim(5165061225, 5165061425)
    ax1.set_xlim(t[294], t[568])
    plt.legend()

    ax1 = fig.add_subplot(313)
    plot, = ax1.plot(t, phi_peak, marker='.', ls="None", color="grey", label="PI")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("RM (rad m$^{-2}$)")
    #ax1.set_xlim(5165061225, 5165061425)
    ax1.set_xlim(t[294], t[568])
    
    fig.savefig("figure.png", bbox_inches="tight")
    plt.close()

# Plot the time-series light curve averaged over the band
def plot_lc(times, I, Q, U, V, nsigma, ppol = "IQUVP"):
    P = np.sqrt(Q*Q + U*U)
    averI = np.nanmean(I, axis=0) * 1000.0  # Average DS in frequency and convert to mJy
    averQ = np.nanmean(Q, axis=0) * 1000.0
    averU = np.nanmean(U, axis=0) * 1000.0
    averV = np.nanmean(V, axis=0) * 1000.0

    averP = np.nanmean(P, axis=0) * 1000.0

    stdI = np.nanstd(averI)
    
    fig = plt.figure(figsize=(7, 5))
    ax1 = fig.add_subplot(111)
    if ppol.find("I") != -1:
        medI = medfilt(averI, 11)
        plot, = ax1.plot(times, averI, marker='', color="black", label="I")
    if ppol.find("Q") != -1:
        medQ = medfilt(averQ, 11)
        plot, = ax1.plot(times, averQ, marker='', color="red", label="Q")
    if ppol.find("U") != -1:
        medU = medfilt(averU, 11)
        plot, = ax1.plot(times, averU, marker='', color="blue", label="U")
    if ppol.find("V") != -1:
        medV = medfilt(averV, 11)
        plot, = ax1.plot(times, averV, marker='', color="green", label="V")
    if ppol.find("P") != -1:
        medP = medfilt(averP, 11)
        plot, = ax1.plot(times, averP, marker='', color="gray", label="P")
    ax1.set_title("Light Curve")
    ax1.set_xlabel("Integration")
    ax1.set_ylabel("Flux Density (mJy)")
    ax1.set_xlim(times[294], times[568])
    ax1.set_ylim(-nsigma*stdI/6.0, nsigma*stdI)
    plt.legend()
    plt.tight_layout()
    fig.savefig("figure2.png", bbox_inches="tight")
#    plt.show()
    plt.close()

def plot_ds(times, freqs, It, Qt, Ut, Vt, sigma, real_time=False, real_freq=False):
    vstd = np.nanstd(Vt)

    current_cmap = cm.get_cmap("cubehelix").copy()
    current_cmap.set_bad(color=current_cmap(0.5))
    nint = len(times)
    nchan = len(freqs)
    if real_time==False and real_freq==False:
        ext = [0, nint, 0, nchan]
        tlabel = "Integration"
        flabel = "Chan"
    elif real_time==True and real_freq==False:
        ext = [0.0, times[-1]-times[0], 0, nchan]
        tlabel = "Elapsed Time (s)"
        flabel = "Chan"
    elif real_time==False and real_freq==True:
        tlabel = "Integration"
        flabel = "$\\nu$ (GHz)"
        ext = [0, nint, freqs[0],freqs[-1]]
    else:
        ext = [0.0, times[-1]-times[0], freqs[0],freqs[-1]]
        tlabel = "Elapsed Time (s)"
        flabel = "$\\nu$ (GHz)"

    fig = plt.figure(figsize=(14, 8))
    ax1 = fig.add_subplot(221)
    ax1.set_title("Dynamic Spectra (I)")
    ax1.set_xlabel(tlabel)
    ax1.set_ylabel(flabel)
    plt.imshow(It, origin='lower', clim=(0.0, sigma*vstd), extent=ext, aspect="auto", cmap=current_cmap)
    cbar = plt.colorbar()
    cbar.set_label('Jy')#,labelpad=-75)

    ax2 = fig.add_subplot(222)
    ax2.set_title("Dynamic Spectra (Q)")
    ax2.set_xlabel(tlabel)
    ax2.set_ylabel(flabel)
    plt.imshow(Qt, origin='lower', clim=(-sigma*vstd, sigma*vstd), extent=ext, aspect="auto", cmap=current_cmap)
    cbar = plt.colorbar()
    cbar.set_label('Jy')#,labelpad=-75)

    ax3 = fig.add_subplot(223)
    ax3.set_title("Dynamic Spectra (U)")
    ax3.set_xlabel(tlabel)
    ax3.set_ylabel(flabel)
    plt.imshow(Ut, origin='lower', clim=(-sigma*vstd, sigma*vstd), extent=ext, aspect="auto", cmap=current_cmap)
    cbar = plt.colorbar()
    cbar.set_label('Jy')#,labelpad=-75)

    ax4 = fig.add_subplot(224)
    ax4.set_title("Dynamic Spectra (V)")
    ax4.set_xlabel(tlabel)
    ax4.set_ylabel(flabel)
    plt.imshow(Vt, origin='lower', clim=(-sigma*vstd, sigma*vstd), extent=ext, aspect="auto", cmap=current_cmap)
    cbar = plt.colorbar()
    cbar.set_label('Jy')#,labelpad=-75)

    plt.tight_layout()
    plt.show()
    plt.close()

def get_fdft(t, freqs, It, Qt, Ut, Vt, startPhi = -1000.0, dPhi = 1.0):
    stopPhi = -startPhi+dPhi

    return getFDF(Qt[t], Ut[t], freqs, startPhi, stopPhi, dPhi)

def plot_sed(t, freqs, i, q, u, v):
    it = i[:,t] * 1000.0  # Average DS in frequency and convert to mJy
    qt = q[:,t] * 1000.0
    ut = u[:,t] * 1000.0
    vt = v[:,t] * 1000.0
    pt = np.sqrt(qt*qt + ut*ut)
    print(np.nanmin(pt), np.nanmax(pt), np.nanmean(pt))

    fig = plt.figure(figsize=(7, 5))
    ax1 = fig.add_subplot(111)
    plots = []
    #plot, = ax1.plot(aver_f, averXX, marker='', color="green", label="XX")
    #plot, = ax1.plot(aver_f, averYY, marker='', color="blue", label="YY")
    plot, = ax1.plot(freqs, it, marker='', color="black", label="I")
    plot, = ax1.plot(freqs, qt, marker='', color="red", label="Q")
    plot, = ax1.plot(freqs, ut, marker='', color="blue", label="U")
    #plot, = ax1.plot(aver_f, modelfit, color="red", lw=2, label="model")
    plot, = ax1.plot(freqs, pt, marker='', color="brown", label="P")
    plot, = ax1.plot(freqs, vt, marker='', color="green", label="V")
    ax1.set_xlabel("Frequency (MHz)")
    ax1.set_ylabel("Flux Density (mJy)")
    plt.legend()
    plt.show()
    
ix = xarray.open_dataset("1658342033_sdp_l0_GPM1839-10_polcal_scan7.ms_StokesI-kata.nc")
qx = xarray.open_dataset("1658342033_sdp_l0_GPM1839-10_polcal_scan7.ms_StokesQ-kata.nc")
ux = xarray.open_dataset("1658342033_sdp_l0_GPM1839-10_polcal_scan7.ms_StokesU-kata.nc")
vx = xarray.open_dataset("1658342033_sdp_l0_GPM1839-10_polcal_scan7.ms_StokesV-kata.nc")

i = ix["I"].values
q = qx["Q"].values
u = ux["U"].values
v = vx["V"].values
t = ix.time.values
print(np.median(t[1:]-t[:-1]))
freqs = ix.freq.values
#plot_ds(t, freqs, i, q, u, v, 10.0, True, True)
#plot_lc(t, i, q, u, v, 20.0)
#plot_sed(418, freqs, i, q, u, v)
#tt=418
#fout = open("spectra_%04d.txt" %(tt), "wt")
#for index in range(len(freqs)):
#    fout.write("%f %f %f %f %f %f %f\n" %(freqs[index], i[index,tt], q[index,tt], u[index,tt], 2.0, 1.5, 1.5))
#fout.close()

plot_fdf_2D(t, freqs, np.transpose(i), np.transpose(q), np.transpose(u), np.transpose(v), -1000.0, 1.0 )
