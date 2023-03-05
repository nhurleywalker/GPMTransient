# coding: utf-8
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

arr = np.loadtxt("t_pa_with_err.txt")
t = arr.T[0]
pa = arr.T[1]

# Loads of gaps so do nearest-neighbour interpolation
dt = t[1] - t[0]
# Need to introduce a step just before this OPM, otherwise it gets interpolated over 10s of seconds
pa[-132] = 50
t = t[pa != 0.0]
pa = pa[pa != 0.0]
f = interp1d(t, pa, kind='nearest')

orig_pa = pa

x = np.arange(t[0], t[-1], dt)
pa = f(x)

# Detect OPMs
lt = np.squeeze(np.where(pa < 0))
steps = lt[1:-1] - lt[0:-2]

# Add dummy entry so that last OPM is detected
steps = np.append(steps, 2)

# Find where the jumps are
ends = np.argwhere(steps > 1)
starts = ends + 1
starts = np.append(0, starts)

# Loop over segments and find location in time and phase
opms = np.empty(len(ends))
opmt = np.empty(len(ends))
for i in range(0, len(ends)):
    a = np.squeeze(lt[starts[i]])
    b = np.squeeze(lt[ends[i]])
    opms[i] = np.mean(pa[a:b])
    opmt[i] = 0.065*np.nanmean([a, b])

cut = 4
diffs = opmt[1:] - opmt[0:-1]
retain = diffs[diffs < cut]
avg_diff = np.nanmean(retain)

fix = np.squeeze(np.argwhere(diffs > cut))

# Don't use the last value because it's not right
for i in fix[0:-1]:
    retain = np.append(retain, diffs[i] / round(diffs[i]/avg_diff))

p = np.nanmean(retain)

pt = opmt[0] + np.arange(0, np.max(0.065 * np.arange(0, len(pa))), p)


fig = plt.figure(figsize=(10, 3))
ax = fig.add_subplot(111)
ax.set_xlabel("time / s")
ax.set_ylabel("PA / deg")
for i in pt:
    ax.axvline(i, color='grey', zorder=-1, alpha=0.5)
ax.scatter(0.065 * np.arange(0, len(pa)), pa, marker=".", color='blue')
ax.scatter(opmt, opms, marker="x", color="red", alpha=0.5)
ax.set_title("Cadence = {0:2.2f}s".format(p))
fig.savefig("pa_wrt_time_htru.png", bbox_inches="tight")

