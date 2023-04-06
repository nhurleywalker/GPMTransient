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
min_nch_frac = 0.25 # Only include time bins with more channels used to calculate them than this fraction of the max number of channels. Example, if a lightcurve is constructed from a total of 100 channels, then only consider time bins that included at least 25 valid channels

table = []

for yaml_file in yaml_files:

    obsid = yaml_file[:10]

    # Get the corresponding TOA
    # (assumes they've already been made before running this script)
    tim_file = f"ppdot_search/{obsid}_orig.tim"
    with open(tim_file, 'r') as tf:
        line = tf.readlines()[0] # These files only have one line in them
    mjd_str = line.split()[2]
    toa = Time(mjd_str, format='mjd')

    with open(yaml_file, 'r') as yf:
        params = parse_yaml(yf) # returns dictionary

    dynspec = Dynspec(**params)

    coord = SkyCoord(ra=params['RA']*u.hr, dec=params['Dec']*u.deg, frame='icrs')
    time = Time(dynspec.t[0] - dynspec.dt/2, format='gps')
    utc = time.utc.datetime.strftime("%Y-%m-%d %H:%M")

    # Get barycentric correction
    bc_correction = TimeDelta(bc_corr(coord, time, EPHEMERIS), format='sec')
    toa += bc_correction

    # Get a few bits of metadata for the table
    telescope = params['telescope']

    midfreq = f"{0.5*(dynspec.f[0] + dynspec.f[-1]):.2f}"

    # Get the corresponding lightcurve files
    # (assumes they've already been made before running this script)
    lightcurve_file = f"{obsid}_lightcurve.txt"
    lightcurve_data = np.loadtxt(lightcurve_file)

    new_pulse_number = round((toa.gps - int(yaml_files[0][:10]))/P + 0.15) + 1 # <--- 0.15 is a rough, "manual" pulse centering

    if new_pulse_number != prev_pulse_number:

        # If this is the first of two (or more) observations in the same pulse, there should be a .tim file
        # for the joint TOA
        try:
            tim_file = f"ppdot_search/{obsid}_mod.tim"
            with open(tim_file, 'r') as tf:
                line = tf.readlines()[0] # These files only have one line in them
            mjd_str = line.split()[2]
            toa = Time(mjd_str, format='mjd')
            bc_correction = TimeDelta(bc_corr(coord, time, EPHEMERIS), format='sec')
            toa += bc_correction
        except:
            pass

        # Parse the light curve
        lc = lightcurve_data[:,1]
        t = lightcurve_data[:,0]
        t_idxs = list(range(len(t)))
        nch = lightcurve_data[:,2] # The number of channels that contributed to each point in the lightcurve
        dt = t[1] - t[0]
        fluence_bins = lc * dt

        # Get the peak flux
        peak_flux_density = np.max(lc)

        # Get the total fluence (only include bins where the relative noise is no more than 2*min,
        # which occurs when number of channels included in lightcurve calc is < 4*max
        fluence = np.sum([fluence_bins[i] for i in range(len(fluence_bins)) if nch[i] >= max(nch)*min_nch_frac])

        row = {
            'num_obs': 1,
            'pulse_number': new_pulse_number,
            'utcs': [utc],
            'toa': f"{toa.mjd:.4f}",
            'telescopes': [telescope],
            'midfreq': [midfreq],
            'peak': f"{peak_flux_density:.0f}",
            'fluence': f"{fluence:.0f}",
        }

    else:
        # Pop the previous row from the table
        row = table.pop()

        # Update row values
        row['num_obs'] += 1
        row['utcs'].append(utc)
        row['telescopes'].append(telescope)
        row['midfreq'].append(midfreq)

        # For the peak and fluence, we have to update the lightcurve first

        # Combine the light curve with the previous one
        # "Combine" means a weighted sum (equivalent to forming the light curve from both dynamic spectra)
        new_t = lightcurve_data[:,0]
        new_lc = lightcurve_data[:,1]
        new_nch = lightcurve_data[:,2] # The number of channels that contributed to each point in the lightcurve
        new_dt = new_t[1] - new_t[0]
        new_fluence_bins = new_lc * new_dt

        for new_i in range(len(new_t)): # Loop through the new time axis
            t_idx = round((new_t[new_i] - t[0])/dt) # The bin number relative to the previous lightcurve
            if t_idx in t_idxs: # If this bin overlaps with the previous lightcurve
                i = t_idxs.index(t_idx)
                lc[i] = (lc[i]*nch[i]*dt + new_lc[new_i]*new_nch[new_i]*new_dt) / (nch[i]*dt + new_nch[new_i]*new_dt) # Weighted sum
                nch[i] += new_nch[new_i] # Update the channel count
                fluence_bins[i] = (fluence_bins[i]*dt + new_fluence_bins[new_i]*new_dt) / (nch[i] + new_nch[new_i])
            else: # If this bin does not overlap with the previous lightcurve, then none of the others will either, so just append the rest in bulk
                t = np.hstack((t, new_t[new_i:]))
                lc = np.hstack((lc, new_lc[new_i:]))
                nch = np.hstack((nch, new_nch[new_i:]))
                fluence_bins = np.hstack((fluence_bins, new_fluence_bins[new_i:]))
                break

        # Same peak and fluence calculations as before
        peak_flux_density = np.max(lc)
        fluence = np.sum([fluence_bins[i] for i in range(len(fluence_bins)) if nch[i] >= max(nch)*min_nch_frac])

        row['peak'] = f"{peak_flux_density:.0f}"
        row['fluence'] = f"{fluence:.0f}"

    # Add the row to the table
    table.append(row)
    prev_pulse_number = new_pulse_number

# Format the table to LaTeX table format
for row in table:
    def multirow(val):
        return f"\multirow{{{row['num_obs']}}}{{*}}{{{val}}}"

    if row['num_obs'] == 1:
        print(f"{row['pulse_number']} & {row['utcs'][0]} & {row['toa']} & {row['telescopes'][0]} & {row['midfreq'][0]} & {row['peak']} & {row['fluence']} \\\\")
    else:
        print(f"{multirow(row['pulse_number'])} & {row['utcs'][0]} & {multirow(row['toa'])} & {row['telescopes'][0]} & {row['midfreq'][0]} & {multirow(row['peak'])} & {multirow(row['fluence'])} \\\\")
        for i in range(1, row['num_obs']):
            print(f"  & {row['utcs'][0]} & {multirow(row['toa'])} & {row['telescopes'][0]} & {row['midfreq'][0]} &   &   \\\\")
