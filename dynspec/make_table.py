'''
Compiles a bunch of information from the yaml files and the light curves,
and prepares it into a table structure for inclusion in the paper.

The output table has the columns:
    1. start time of obs in UTC
    2. telescope
    3. mid-frequency
    4. Barycentred TOA -- dedispersed to infinite frequency
    5. brightest flux density recorded (peak pixel)
    6. fluence (sigma-clipped at 3 sigma)
    7. corrected flux density (e.g. extrapolated to 1 GHz)

Before running this script, make sure that the lightcurves are available:
    make all_lightcurves

Also make sure the TOA analysis has already been done:
    (See the README in the ppsot_search folder)
'''

import yaml
import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.stats import sigma_clip
import datetime
from dedisperse_dynspec import Dynspec, parse_yaml
from bc_corr import bc_corr

EPHEMERIS = 'de430.bsp'

yaml_files = glob.glob('*.yaml')

table = []

for yaml_file in yaml_files:

    # Get the corresponding lightcurve files
    # (assumes they've already been made before running this script)
    obsid = yaml_file[:10]
    lightcurve_file = f"{obsid}_lightcurve.txt"
    lightcurve_data = np.loadtxt(lightcurve_file)
    lc = lightcurve_data[:,1]
    t = lightcurve_data[:,0]
    dt = t[1] - t[0]

    # Now we just want the bits with signal in it, so do a first round of sigma clipping
    result = sigma_clip(lc, sigma=3, masked=True, sigma_lower=np.inf)
    signal = lc[result.mask]
    clipped_t = t[result.mask]

    if len(clipped_t) > 0: # If SOME signal was found...

        # Define "real" signal to be any time bin within some interval (e.g. 10 sec)
        # from a time bin that survived the sigma clipping.
        interval = 10
        expanded_idxs = [idx for idx in range(len(t)) if np.min(np.abs(t[idx] - clipped_t)) <= interval]
        expanded_signal = lc[expanded_idxs]
        expanded_t = t[expanded_idxs]
    else: # Maybe the thing is too noisy. For now, just include everything
        expanded_t = t
        expanded_signal = lc

    plt.clf()
    ax = plt.gca()
    plt.fill_between(t, 0, 1, where=[t0 in expanded_t for t0 in t], color='gray', alpha=0.5, transform=ax.get_xaxis_transform())
    plt.plot(t, lc)
    #plt.plot(expanded_t, expanded_signal)
    plt.xlabel("Time (s)")
    plt.ylabel("Flux density (Jy)")
    plt.savefig(f"{obsid}_sigmaclip.png")

    # Get the corresponding TOA
    tim_file = f"ppdot_search/{obsid}.tim"
    with open(tim_file, 'r') as tf:
        line = tf.readlines()[0] # These files only have one line in them
        mjd_str = line.split()[2]
        toa = Time(mjd_str, format='mjd')

    # Get the peak flux
    peak_flux_density = np.max(lightcurve_data[:,1])

    with open(yaml_file, 'r') as yf:
        params = parse_yaml(yf) # returns dictionary
        dynspec = Dynspec(**params)

        # Get barycentric correction
        coord = SkyCoord(ra=params['RA']*u.hr, dec=params['Dec']*u.deg, frame='icrs')
        time = Time(dynspec.t[0] - dynspec.dt/2, format='gps')
        bc_correction = bc_corr(coord, time, EPHEMERIS)

        row = {
            "utc": time.utc.datetime.strftime("%Y-%m-%d %H:%M"),
            "telescope": params['telescope'],
            "midfreq": f"{0.5*(dynspec.f[0] + dynspec.f[-1]):.2f}",
            "peak": peak_flux_density,
            "toa": toa.mjd + bc_correction/86400,
        }
        table.append(row)

#print(table)
