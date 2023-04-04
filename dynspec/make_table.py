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

'''

import yaml
import numpy as np
import glob
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
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
            "toa": toa.mjd,
        }
        table.append(row)

print(table)
