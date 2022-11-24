#!/usr/bin/env python

import os
import numpy as np
from glob import glob
from matplotlib import pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf
from astropy.time import Time
import pandas as pd
import yaml

import matplotlib.font_manager
from matplotlib import rc
# Nature requires sans-serif fonts
plt.rcParams.update({
    "text.usetex": False,
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

# inches to cm
cm = 1/2.54

# Just the auto-correlations for selected pulses, for the paper

# https://stackoverflow.com/questions/67895987/how-to-add-two-partially-overlapping-numpy-arrays-and-extend-the-non-overlapping
def avg_signal(a, b, ai=0, bi=0):
    assert ai >= 0
    assert bi >= 0

    al = len(a)
    bl = len(b)
    cl = max(ai + al, bi + bl)
    c = np.zeros(cl)
    c[ai: ai + al] += a
    c[bi: bi + bl] += b
    # unweighted average
    c[bi: ai + al] /= 2
    return c

# https://stackoverflow.com/questions/24838629/round-off-float-to-nearest-0-5-in-python
def rr(number):
    """Round a number to the closest half integer.
    >>> round_off_rating(1.3)
    1.5
    >>> round_off_rating(2.6)
    2.5
    >>> round_off_rating(3.0)
    3.0
    >>> round_off_rating(4.1)
    4.0"""

    return np.round(number * 2) / 2

def sc(data, level=3, loops=4):
    ''' some sigma-clipping '''
    std = np.nanstd(data)
    for i in range(0, loops):
        std = np.nanstd(data[data<level*std])
    return std

def make_sparse(array, start=0, box=5):
    ''' remove close-packed points '''
    i = 0
    keep = []
    print(start)
    for j in range(0, len(array)-1):
        if array[j] > start:
            diff = array[j+1] - array[j]
            if diff < box:
                i += 1
            else:
                keep.append(array[j - int(i/2) - 1])
                i = 0
    return(keep)
        

P = rr(1318.19578) #Period in seconds
tscat200 = 0.5 # Scattering at 200MHz in seconds

selection = ["1341931472", "1342379534", "1342623496"]

xlims = [[20, 250], [90, 350], [100, 280]]
tstarts = [20, 90, 100]
nsec = [230, 260, 180]
# Bit hand-holdy but the computer isn't great at this
acfstarts = [40, 0, 50]
cutoffs = [-0.29, -0.20, -0.26]


cs = []
path = "../dynspec/ppdot_search/"
for s in selection:
    cs.append(path+s+"_lightcurve.txt")

print("Analysing {0} light curves".format(len(cs)))

# https://www.nature.com/nature/for-authors/final-submission#:~:text=For%20guidance%2C%20Nature's%20standard%20figure,(120%E2%80%93136%20mm).
fig = plt.figure(figsize=(17.8*cm,10*cm))

for i in range(0, len(selection)):
    arr1 = np.loadtxt(cs[i])
    obsid1 = os.path.split(cs[i])[-1][0:10]
    ti = Time(obsid1, format="gps")
    ti.format="ymdhms"
    td = ti.value

    yaml_file = "../dynspec/{0}.yaml".format(obsid1)
    with open(yaml_file, "r") as stream:
        yaml_params = yaml.safe_load(stream)
    nu = float(yaml_params['Dynamic spectrum']['Centre of lowest channel (MHz)'])
    ts = float(yaml_params['Dynamic spectrum']['Sample time (s)'])

    s1 = arr1[:,1]
    ns1 = s1 - np.mean(s1)
    var = np.var(s1)
    t1 = arr1[:,0]
    acorr = np.correlate(ns1, ns1, 'full')[len(ns1)-1:]

# Find peaks
    grad = np.gradient(np.gradient(acorr))
    grad = grad/(np.nanmax(grad)- np.nanmin(grad))
# Only find one
    peaks = np.squeeze(np.where(grad < cutoffs[i]))
    print(peaks)
    print(make_sparse(np.append(peaks, 999999), acfstarts[i]))
    peaks = make_sparse(np.append(peaks, 999999), acfstarts[i])
 
    lc_sub = i + 1 
    acf_sub = i + 4

    t = ts*np.arange(0,len(acorr),1)
    ax_lc = fig.add_subplot(2,3,lc_sub)
    ax_lc.plot(t, s1/np.nanmax(s1), alpha=1, lw=0.5, color="blue")
    #ax_lc.set_title(f"{obsid1}: light curve")
    ax_lc.set_title(f"{td[0]:04d}-{td[1]:02d}-{td[2]:02d} {td[3]:02d}:{td[4]:02d}:{td[5]:02.0f}")
    ax_lc.set_xlim([tstarts[i], tstarts[i] + nsec[i]])
    ax_acf = fig.add_subplot(2,3,acf_sub)
    ax_acf.plot(t, acorr/np.nanmax(acorr), alpha=1, lw=0.5, color="blue")
# Helper info to show how the gradient detector works
#    ax_acf.plot(t, grad, alpha=1, lw=0.5, color="red")
#    ax_acf.axhline(cutoffs[i], lw=0.25)
    ax_acf.set_xlim([0, nsec[i]/2])
    try:
        for peak in peaks:
            ax_acf.axvline(t[peak], lw=0.25, color="grey")
            ax_acf.text(t[peak], 0.9, f" {t[peak]:3.0f}s")
    except TypeError:
        ax_acf.axvline(t[peaks], lw=0.25, color="grey")
        ax_acf.text(t[peaks], 0.9, f" {t[peaks]:3.0f}s")
    #ax_acf.set_title(f"{obsid1}: autocorrelation")
    if i < 1:
        ax_lc.set_ylabel("light curve (a.u.)")
        ax_acf.set_ylabel("correlated power (a.u.)")
    else:
        ax_lc.yaxis.set_ticklabels([])
        #ax_acf.yaxis.set_ticklabels([])
    if i == 1:
        ax_lc.set_xlabel("time / seconds")
        ax_acf.set_xlabel("lag / seconds")
helper = "acf.pdf"
fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
fig.savefig(helper, bbox_inches="tight")
#    ax.axvline(50., alpha=0.3, color='k')
#    ax.axvline(100., alpha=0.3, color='k')

    # Auto-correlations
#    for i in range(0, len(final_list)):
#        print(final_cs[i])
#        print(os.path.split(final_cs[i])[-1])
#    obsid1 = os.path.split(final_cs[i])[-1][0:10]
#    yaml_file = "../dynspec/{0}.yaml".format(obsid1)
#    with open(yaml_file, "r") as stream:
#        yaml_params = yaml.safe_load(stream)
#    nu = float(yaml_params['Dynamic spectrum']['Centre of lowest channel (MHz)'])
#    ts = float(yaml_params['Dynamic spectrum']['Sample time (s)'])
#    tscat = tscat200 * (nu / 200.)**-4
#    arr1 = final_list[i]
#    s1 = arr1[:,1]
#    ns1 = s1 - np.mean(s1)
#    var = np.var(s1)
#    t1 = arr1[:,0]
#    acorr = np.correlate(ns1, ns1, 'full')[len(ns1)-1:]
#    acorr = acorr / var / len(ns1)
#    fig = plt.figure()
#    ax = fig.add_subplot(211)
#    ax.plot(ts*np.arange(0,len(acorr),1), s1, alpha=1, lw=1, color="blue")
#    ax.set_ylabel("light curve (a.u.)")
#    ax.set_title(f"{obsid1}: autocorrelation")
#    ax = fig.add_subplot(212)
#    ax.plot(ts*np.arange(0,len(acorr),1), acorr, alpha=1, lw=1, color="blue")
#    ax.set_xlabel("lag / seconds")
#    ax.set_ylabel("correlated power (a.u.)")
#    ax.axvline(50., alpha=0.3, color='k')
#    ax.axvline(100., alpha=0.3, color='k')
#    ax.axvspan(0, tscat, alpha=0.3, color='r')
#    #ax.set_ylim([-0.5, 1.2])
#    helper = f"{obsid1}_acf.png"
#    fig.savefig(helper, bbox_inches="tight")
#    plt.close(fig)

# TODO: highlight the 2Ps in some way

            # Some QA to add back in at some point
            # Some QA to add back in at some point
    #        if np.isnan(s1).any() or np.isnan(s2).any():
    #            print(f"found some nans for {cs[i]} or cs[i+1], not proceeding further")
        #        elif cs[i] == "1342096266" or cs[i+1] == "1342096266":
        #            print("Skipping Parkes data")
    #        else:
