# coding: utf-8
import numpy as np
from matplotlib import pyplot as plt

arr = np.loadtxt("pa_with_err.txt")
pa = arr.T[0]


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

cut = 3
diffs = opmt[1:] - opmt[0:-1]
#print(diffs)
retain = diffs[diffs < cut]
#print(retain)
avg_diff = np.nanmean(retain)

fix = np.squeeze(np.argwhere(diffs > cut))
#print(fix)

for i in fix:
    retain = np.append(retain, diffs[i] / round(diffs[i]/avg_diff))

p = np.nanmean(retain)

pt = opmt[0] + np.arange(0, np.max(0.065 * np.arange(0, len(pa))), p)

fig = plt.figure(figsize=(10, 3))
ax = fig.add_subplot(111)
ax.set_xlabel("time / s")
ax.set_ylabel("PA / deg")
for i in pt:
    ax.axvline(i, color='grey', zorder=-1, alpha=0.5)
ax.scatter(0.065 * np.arange(0, len(pa)), pa, marker=".")
ax.scatter(opmt, opms, marker="x", color="red", alpha=0.5)
ax.set_title("Cadence = {0:2.2f}s".format(p))
fig.savefig("pa_wrt_time_htru.png", bbox_inches="tight")

