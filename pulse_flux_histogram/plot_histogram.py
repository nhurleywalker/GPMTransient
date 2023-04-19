#from astropy.table import Table
#from astropy.coordinates import SkyCoord
#from astropy import unit as u
#import scipy.ndimage
#from scipy.ndimage.filters import gaussian_filter
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.font_manager
from matplotlib import rc
from matplotlib.patches import Rectangle

import pandas as pd

# Nature requires sans-serif fonts
plt.rcParams.update({
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

cm = 1/2.54  # centimeters in inches

df = pd.read_csv("pulse_table.csv")
# Pulse number,UTC,Barycentred TOA (MJD),Telescope,Frequency (MHz),Peak flux density at freq (Jy),Peak flux density at 1 GHz (mJy),Fluence at freq (Jy s),Fluence at 1 GHz (Jy s)

#https://www.nature.com/nature/for-authors/final-submission#:~:text=For%20guidance%2C%20Nature's%20standard%20figure,(120%E2%80%93136%20mm).
fig = plt.figure(figsize=(8.9*cm,7.8*cm))
ax1 = fig.add_subplot(111)
h1 = ax1.hist(df["Peak flux density at 1 GHz (mJy)"][df["Fluence at 1 GHz (Jy s)"]>0], color='red', alpha=0.5, label="$S$")
ax2 = ax1.twiny()
h2 = ax2.hist(df["Fluence at 1 GHz (Jy s)"][df["Fluence at 1 GHz (Jy s)"]>0], color='blue', alpha=0.5, label="Fluence")
ax2.set_xlabel("Fluence$_\mathrm{1GHz}$ (Jy s)")
ax1.set_xlabel("$S_\mathrm{1GHz}$ (mJy)")
ax1.set_ylabel("Number of pulses")
ax1.set_yticks([0, 5, 10, 15, 20])
ax1.yaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
# https://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend
colors = ['red', 'blue']
handles = [Rectangle((0, 0), 1, 1, color=c, alpha=0.5) for c in colors]
labels = ["Fluence", "$S$"]
plt.legend(handles, labels)
fig.savefig("S_hist.pdf", bbox_inches="tight", dpi=300)
