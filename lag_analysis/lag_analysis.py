#!/usr/bin/env python

import os
import numpy as np
from glob import glob
from matplotlib import pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf
import pandas as pd
import yaml

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

P = rr(1318.19578) #Period in seconds
tscat200 = 0.5 # Scattering at 200MHz in seconds

path = "../dynspec/"
cs = sorted(glob(path+"*_lightcurve.txt"))
# Remove Parkes data, since it is on a different time scale and is difficult to average with the MWA data
try:
    cs.remove(path+"1342096266_lightcurve.txt")
except ValueError:
    pass

final_list = []
final_cs = []
# deal with light curves that overlap
skipnext = False
for i in range(0, len(cs)):
    outfile = cs[i].replace(".txt", "_mod.txt")
    if i < len(cs) -1:
        arr1 = np.loadtxt(cs[i])
        t1 = rr(arr1.T[0])
        s1 = arr1.T[1]
        arr2 = np.loadtxt(cs[i+1])
        t2 = rr(arr2.T[0])
        s2 = arr2.T[1]
        obsid1 = os.path.split(cs[i])[-1][0:10]
        yaml_file = "../dynspec/{0}.yaml".format(obsid1)
        with open(yaml_file, "r") as stream:
            yaml_params = yaml.safe_load(stream)
        nu = float(yaml_params['Dynamic spectrum']['Centre of lowest channel (MHz)'])
        ts = float(yaml_params['Dynamic spectrum']['Sample time (s)'])
        if skipnext is False:
        # light curves overlap
            if np.isin(t1,t2).any():
                print(f"merging {cs[i]} and {cs[i+1]}")
        # don't normalise
        # TODO: normalise by the differences in the radio SED
#                rescale = np.nanmean(s1[np.isin(t1,t2)]/s2[np.isin(t2,t1)])
# This gave bad results at times
                tdiff = np.squeeze(np.where(t1 == t1[np.isin(t1,t2)][0]))
                s3 = avg_signal(s1, s2, 0, tdiff)
                #s3 = avg_signal(s1/np.nanmax(s1), s2/np.nanmax(s2), 0, tdiff)
                t3 = avg_signal(t1, t2, 0, tdiff)
                skipnext = True
#                fig = plt.figure()
#                ax = fig.add_subplot(111)
#                ax.plot(t1, s1/np.nanmax(s1), alpha=0.5, lw=2, color="blue", label="t1")
#                ax.plot(t2, s2/np.nanmax(s2), alpha=0.5, lw=2, color="green", label="t2")
#                ax.plot(t3, s3/np.nanmax(s3), alpha=0.5, lw=2, color="red", label="t3")
#                ax.legend()
#                ax.set_title(tdiff)
#                helper = f"merge_{cs[i][0:10]}_{cs[i+1][0:10]}.png"
#                fig.savefig(helper)
        # light curves have a very short gap (<3min)
            elif t2[0] - t1[-1] <= 360:
                print(f"joining {cs[i]} and {cs[i+1]} ({i}th check)")
                #s3 = np.concatenate((s1/np.nanmax(s1), np.zeros(int((t2[0]-t1[-1])/0.5)-1), s2/np.nanmax(s2)))

                pad = np.arange(t1[-1] + ts, t2[0], ts)
                t3 = np.concatenate((t1, pad, t2))
 # Try interpolating
                s3 = np.concatenate((s1, np.array(((s2[0] - s1[-1])/(t2[0] - t1[-1]))*(pad - t1[-1]) + s1[-1]), s2))
# np.zeros(len(pad)), s2))
                skipnext = True
        # unique light curve that stands alone
            else:
                print(f"solo: {cs[i]} ({cs[i+1]} is not connected)")
                #s3 = s1/np.nanmax(s1)
                s3 = s1/np.nanmax(s1)
                t3 = t1
            #combine = np.stack([t3,s3], axis=1)
            # New change: normalise here
            combine = np.stack([t3,s3/np.nanmax(s3)], axis=1)
            final_list.append(combine)
            final_cs.append(cs[i])
            np.savetxt(outfile, combine)
        else:
            skipnext = False

lags = np.arange(-120, 120, 1)
nP = np.zeros(len(final_list))

corr_matrix = np.zeros((len(final_list), len(lags)))

# work out which light curves are one period apart -- that is a pair
for i in range(0, len(final_list)):
    if i < len(final_list) -1:
        arr1 = final_list[i]
        arr2 = final_list[i+1]
        t1 = np.copy(arr1[:,0])
        t2 = np.copy(arr2[:,0])
        s1 = arr1[:,1]
        s2 = arr2[:,1]
        tdiff = t2[0] - t1[-1]
        obsid1 = os.path.split(final_cs[i])[-1][0:10]
        obsid2 = os.path.split(final_cs[i+1])[-1][0:10]
#        print(f"Checking {final_cs[i]} (spans {t1[0]}--{t1[-1]}) and {final_cs[i+1]} (spans {t2[0]}--{t2[-1]}); tdiff = {tdiff}")
        if tdiff < 2*P and tdiff > 0 :
        # Normalise so the time axes are the same
        # This will allow us to crop the data in each lag
            if tdiff < P:
                print(f"{obsid1} and {obsid2} match (tdiff = {tdiff}), subtracting one period")
                t2 -= P
                nP[i] = 1
            else:
                print(f"{obsid1} and {obsid2} nearly match (tdiff = {tdiff}), subtracting two periods")
                t2 -= 2*P
                nP[i] = 2
            for j in range(0, len(lags)):
                t3 = t2 + lags[j]
                # Second data set happens last
                if t3[0] > t1[0] and t3[-1] > t1[-1]:
                    t1_start = np.squeeze(np.where(t1 == t3[0]))
                    t1_end = -1
                    t3_start = 0
                    t3_end = np.squeeze(np.where(t3 == t1[-1]))
                # first data set happens last
                elif t1[0] > t3[0] and t1[-1] > t3[-1]:
                    t1_start = 0
                    t1_end = np.squeeze(np.where(t1 == t3[-1]))
                    t3_start = np.squeeze(np.where(t3 == t1[0]))
                    t3_end = -1
                # first dataset is contained within second dataset
                elif t1[0] > t3[0] and t1[-1] < t3[-1]:
                    t1_start = 0
                    t1_end = -1
                    t3_start = np.squeeze(np.where(t3 == t1[0]))
                    t3_end = np.squeeze(np.where(t3 == t1[-1]))
                # second dataset is contained within first dataset
                else:
                    t1_start = np.squeeze(np.where(t1 == t3[0]))
                    t1_end = np.squeeze(np.where(t1 == t3[-1]))
                    t3_start = 0
                    t3_end = -1
                if (t1_start or t1_start == 0) and (t1_end or t1_end == 0) and (t3_start or t3_start == 0) and (t3_end or t3_end == 0):
                    s1_crop = s1[t1_start:t1_end] #/np.nanmax(s1)
                    s3_crop = s2[t3_start:t3_end] #/np.nanmax(s2)
                    if len(s1_crop) == len(s3_crop):
                        corr_matrix[i, j] = np.correlate(s1_crop, s3_crop)
                    else:
                        print(f"{obsid1} and {obsid2} have different lengths {len(s1_crop)} {len(s3_crop)} for lag {lags[j]}")
                else:
                    print(f"{obsid1} and {obsid2} have no overlapping time for lag {lags[j]}?")
#                fig = plt.figure()
#                ax = fig.add_subplot(111)
#                ax.plot(t1, s1/np.nanmax(s1), alpha=0.5, lw=2, color="blue")
#                ax.plot(t3, s2/np.nanmax(s2), alpha=0.5, lw=2, color="red")
#                try:
#                    ax.plot(t1[t1_start:t1_end], s1_crop, lw=1, color="blue")
#                    ax.plot(t3[t3_start:t3_end], s3_crop, lw=1, color="red")
#                except:
#                    pass
#                helper = f"{final_cs[i][0:10]}_{final_cs[i+1][0:10]}_lag{lags[j]}.png"
#                fig.savefig(helper)
#        else:
#            print(f"{final_cs[i]} and {final_cs[i+1]} do not match (tdiff = {tdiff})")

weights = np.zeros(len(final_list))
best_lags = np.argmax(corr_matrix, axis=1)
for i in range(0, len(final_list)):
    if best_lags[i] != 0:
        arr1 = final_list[i]
        arr2 = final_list[i+1]
        t1 = np.copy(arr1[:,0])
        t2 = np.copy(arr2[:,0])
        s1 = arr1[:,1]
        s2 = arr2[:,1]
        t3 = t2 + lags[best_lags[i]] - nP[i]*P
        weights[i] = 1/np.sqrt(sc(s1)**2 + sc(s2)**2)
        obsid1 = os.path.split(final_cs[i])[-1][0:10]
        obsid2 = os.path.split(final_cs[i+1])[-1][0:10]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(t1, s1, alpha=0.5, lw=1, color="blue", label=obsid1)
        ax.plot(t3, s2, alpha=0.5, lw=1, color="red", label=obsid2)
        ax.legend(loc="upper right")
        ax.set_ylim([-0.5, 1.2])
        ax.set_title(f"pair {i}: Optimum lag = {lags[best_lags[i]]} s")
        helper = f"{obsid1}_{obsid2}_optimum_lag{lags[best_lags[i]]}.png"
        fig.savefig(helper, bbox_inches="tight")

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.imshow(np.tile(weights, (corr_matrix.shape[1],1)).T *corr_matrix, origin="lower", extent = [lags[0], lags[-1], 0, len(final_list)], aspect="auto", interpolation="none")
ax.set_ylabel("observation pair")
ax.set_xlabel("lag / seconds")
fig.savefig("correlation_matrix.png", bbox_inches="tight")

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.plot(lags, np.average(corr_matrix, axis=0, weights = weights))
ax.set_xlabel("lag / seconds")
ax.set_ylabel("summed correlation coefficients (a.u.) ")
fig.savefig("correlation_matrix_average.png", bbox_inches="tight")

final_best = lags[np.argmax(np.average(corr_matrix, axis=0, weights = weights))]
print(final_best)

# Try plotting all the light curves with the globally optimised lag
for i in range(0, len(final_list)):
    if best_lags[i] != 0:
        obsid1 = os.path.split(final_cs[i])[-1][0:10]
        obsid2 = os.path.split(final_cs[i+1])[-1][0:10]
        arr1 = final_list[i]
        arr2 = final_list[i+1]
        t1 = np.copy(arr1[:,0])
        t2 = np.copy(arr2[:,0])
        s1 = arr1[:,1]
        s2 = arr2[:,1]
        t3 = t2 + nP[i]*final_best - nP[i]*P
        weights[i] = 1/np.sqrt(sc(s1)**2 + sc(s2)**2)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(t1, s1, alpha=0.5, lw=1, color="blue", label=f"{obsid1}")
        ax.plot(t3, s2, alpha=0.5, lw=1, color="red", label=f"{obsid2}")
        ax.legend(loc="upper right")
        ax.set_ylim([-0.5, 1.2])
        ax.set_title(f"pair {i}: Applied global optimum {final_best} s")
        helper = f"{obsid1}_{obsid2}_force_{final_best}.png"
        fig.savefig(helper, bbox_inches="tight")

fig = plt.figure(figsize=(6,6))

# Auto-correlations
for i in range(0, len(final_list)):
    obsid1 = os.path.split(final_cs[i])[-1][0:10]
    yaml_file = "../dynspec/{0}.yaml".format(obsid1)
    with open(yaml_file, "r") as stream:
        yaml_params = yaml.safe_load(stream)
    nu = float(yaml_params['Dynamic spectrum']['Centre of lowest channel (MHz)'])
    ts = float(yaml_params['Dynamic spectrum']['Sample time (s)'])
    tscat = tscat200 * (nu / 200.)**-4
    arr1 = final_list[i]
    s1 = arr1[:,1]
    ns1 = s1 - np.mean(s1)
    var = np.var(s1)
    t1 = arr1[:,0]
    acorr = np.correlate(ns1, ns1, 'full')[len(ns1)-1:]
    acorr = acorr / var / len(ns1)
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(ts*np.arange(0,len(acorr),1), s1, alpha=1, lw=1, color="blue")
    ax.set_ylabel("light curve (a.u.)")
    ax.set_title(f"{obsid1}: autocorrelation")
    ax = fig.add_subplot(212)
    ax.plot(ts*np.arange(0,len(acorr),1), acorr, alpha=1, lw=1, color="blue")
    ax.set_xlabel("lag / seconds")
    ax.set_ylabel("correlated power (a.u.)")
    ax.axvline(50., alpha=0.3, color='k')
    ax.axvline(100., alpha=0.3, color='k')
    ax.axvspan(0, tscat, alpha=0.3, color='r')
    #ax.set_ylim([-0.5, 1.2])
    helper = f"{obsid1}_acf.png"
    fig.savefig(helper, bbox_inches="tight")
    plt.close(fig)

# TODO: highlight the 2Ps in some way

            # Some QA to add back in at some point
            # Some QA to add back in at some point
    #        if np.isnan(s1).any() or np.isnan(s2).any():
    #            print(f"found some nans for {cs[i]} or cs[i+1], not proceeding further")
        #        elif cs[i] == "1342096266" or cs[i+1] == "1342096266":
        #            print("Skipping Parkes data")
    #        else:
