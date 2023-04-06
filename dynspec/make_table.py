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
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.stats import sigma_clip
import datetime
from dedisperse_dynspec import Dynspec, parse_yaml
from bc_corr import bc_corr

EPHEMERIS = 'de430.bsp'

yaml_files = glob.glob('*.yaml')
yaml_files.sort()

P = 1318.1956 # Approximate period in seconds
prev_pulse_number = None # For keeping track if a given lightcurve is in the same pulse as the previous light curve

table = []

for yaml_file in yaml_files:

    obsid = yaml_file[:10]

    # Get the corresponding TOA
    # (assumes they've already been made before running this script)
    tim_file = f"ppdot_search/{obsid}.tim"
    with open(tim_file, 'r') as tf:
        line = tf.readlines()[0] # These files only have one line in them
        mjd_str = line.split()[2]
        toa = Time(mjd_str, format='mjd')

    new_pulse_number = round((toa.gps - int(yaml_files[0][:10]))/P + 0.15) + 1, # <--- 0.15 is a rough, "manual" pulse centering

    if new_pulse_number != prev_pulse_number:

        # Reset the fluence calc
        fluence = 0

        # Get the corresponding lightcurve files
        # (assumes they've already been made before running this script)
        lightcurve_file = f"{obsid}_lightcurve.txt"
        lightcurve_data = np.loadtxt(lightcurve_file)
        lc = lightcurve_data[:,1]
        t = lightcurve_data[:,0]
        nch = lightcurve_data[:,2] # The number of channels that contributed to each point in the lightcurve
        dt = t[1] - t[0]

        # Get the peak flux
        peak_flux_density = np.max(lightcurve_data[:,1])

        with open(yaml_file, 'r') as yf:
            params = parse_yaml(yf) # returns dictionary

        dynspec = Dynspec(**params)

        # Get barycentric correction
        coord = SkyCoord(ra=params['RA']*u.hr, dec=params['Dec']*u.deg, frame='icrs')
        time = Time(dynspec.t[0] - dynspec.dt/2, format='gps')
        bc_correction = TimeDelta(bc_corr(coord, time, EPHEMERIS), format='sec')
        toa += bc_correction

        row = {
            "utc": time.utc.datetime.strftime("%Y-%m-%d %H:%M"),
            "toa": toa.mjd,
            "pulse_number": new_pulse_number,
            "telescope": params['telescope'],
            "midfreq": f"{0.5*(dynspec.f[0] + dynspec.f[-1]):.2f}",
            "peak": peak_flux_density,
        }
    else:
        table.append(row)

    print(row)
