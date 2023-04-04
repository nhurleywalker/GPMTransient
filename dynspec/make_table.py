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
'''

import yaml
import numpy as np
import glob
from astropy.time import Time
import datetime
from dedisperse_dynspec import Dynspec, parse_yaml

yaml_files = glob.glob('*.yaml')

table = []

for yaml_file in yaml_files:

    # Get the corresponding lightcurve files
    # (assumes they've already been made before running this script)
    lightcurve_file = f"{yaml_file[:10]}_lightcurve.txt"
    lightcurve_data = np.loadtxt(lightcurve_file)

    # Get the peak flux
    #...

    with open(yaml_file, 'r') as yf:
        params = parse_yaml(yf) # returns dictionary
        dynspec = Dynspec(**params)

        row = {
            "utc": Time(dynspec.t[0] - dynspec.dt/2, format='gps').utc.datetime.strftime("%Y-%m-%d %H:%M"),
            "telescope": params['telescope'],
            "midfreq": f"{0.5*(dynspec.f[0] + dynspec.f[-1]):.2f}",
        }
        table.append(row)

print(table)
