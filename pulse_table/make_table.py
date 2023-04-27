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
    (See the README in the ppdot_search folder)
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
import sys
sys.path.append("../dynspec")
from dedisperse_dynspec import Dynspec, parse_yaml
from bc_corr import bc_corr

def round_to_sf(x, sf=3):
    i = int(np.floor(np.log10(abs(x))))
    rounded = round(x, sf - 1 - i)
    if i - sf >= -1:
        return f"{rounded:.0f}"
    else:
        return f"{rounded}"

sf = 3 # Number of significant figures in the table

def curved_law(nu, ref_nu=1*u.GHz, s_nu=149*u.mJy, alpha=-3.17, q=-0.56):
    return s_nu * (nu/ref_nu) ** alpha * \
            np.exp(q * np.log(nu/ref_nu)**2)


EPHEMERIS = '../dynspec/de430.bsp'
TOAS_BARY_PATH = "../toas_bary"
TOAS_NOBARY_PATH = "../toas_nobary"
LIGHTCURVES_BARY_PATH = "../lightcurves_bary"

SHOW_PLOTS = False

yaml_files = glob.glob(f'{LIGHTCURVES_BARY_PATH}/*.yaml')
#yaml_files = [f'{LIGHTCURVES_BARY_PATH}/1342623496.yaml', f'{LIGHTCURVES_BARY_PATH}/1342623792.yaml']  # <-- Just the MWA "split" pulse

yaml_files.sort()

# A hack to get the Murriyang yaml out from the middle of the two MWA ones of the same pulse
a, b = yaml_files.index(f'{LIGHTCURVES_BARY_PATH}/1342096266.yaml'), yaml_files.index(f'{LIGHTCURVES_BARY_PATH}/1342096400.yaml')
yaml_files[b], yaml_files[a] = yaml_files[a], yaml_files[b]

P = 1318.1956 # Approximate period in seconds
prev_pulse_number = None # For keeping track if a given lightcurve is in the same pulse as the previous light curve
min_bw_frac = 0.25 # Only include time bins with more channels used to calculate them than this fraction of total bandwidth. For example, if a lightcurve is constructed from 100 MHz, and if this number is 0.25, then only consider time bins that included at least 25 MHz of bandwidth

table = []
fluence_bins_to_plot = []

for yaml_file in yaml_files:

    print(f"Reading {yaml_file}...")
    obsid = yaml_file[-15:-5]

    # Get the corresponding TOA
    # (assumes they've already been made before running this script)
    tim_file = f"{TOAS_BARY_PATH}/{obsid}.tim"
    with open(tim_file, 'r') as tf:
        line = tf.readlines()[0] # These files only have one line in them
    mjd_str = line.split()[2]
    toa = Time(mjd_str, format='mjd')

    with open(yaml_file, 'r') as yf:
        params = parse_yaml(yf) # returns dictionary

    dynspec = Dynspec(**params)

    coord = SkyCoord(ra=params['RA']*u.hr, dec=params['Dec']*u.deg, frame='icrs')
    time = Time(dynspec.t[0] - dynspec.dt/2, format='gps')
    
    # Old way, which truncates the minute down
    #utc = Time(int(obsid), format='gps').utc.datetime.strftime("%Y-%m-%d %H:%M")
    # Copying what is done in make_pulsestack, which rounds to the nearest minute
    utc_time = Time(int(obsid), format="gps")
    utc_time.format="ymdhms"
    utc_value = utc_time.value
    utc = f"{utc_value[0]:04d}-{utc_value[1]:02d}-{utc_value[2]:02d} {utc_value[3]:02d}:{utc_value[4]:02d}"

    # Get barycentric correction
    #bc_correction = TimeDelta(bc_corr(coord, time, EPHEMERIS), format='sec')
    #toa += bc_correction

    # Get a few bits of metadata for the table
    telescope = params['telescope']

    midfreq = 0.5*(dynspec.f[0] + dynspec.f[-1])

    # Get the corresponding lightcurve files
    # (assumes they've already been made before running this script)
    lightcurve_file = f"{LIGHTCURVES_BARY_PATH}/{obsid}_lightcurve.txt"
    lightcurve_data = np.loadtxt(lightcurve_file)

    new_pulse_number = round((toa.gps - int(yaml_files[0][-15:-5]))/P + 0.15) + 1 # <--- 0.15 is a rough, "manual" pulse centering

    if new_pulse_number != prev_pulse_number or telescope == "Murriyang":

        # If this is the first of two (or more) observations in the same pulse, there should be a .tim file
        # for the joint TOA
        try:
            tim_file = f"{TOAS_NOBARY_PATH}/{obsid}_mod.tim"
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
        bw = params['bw']
        freq_range = nch*bw
        fluence_bins = lc * dt

        # Get the scaling factor for 1 GHz
        scale_factor = (curved_law(1*u.GHz) / curved_law(midfreq*u.MHz)).decompose()

        # Get the peak flux
        peak_flux_density = np.max(lc)

        # Get the total fluence (only include bins where the relative noise is no more than 2*min,
        # which occurs when number of channels included in lightcurve calc is < 4*max
        fluence = np.sum([fluence_bins[i] for i in range(len(fluence_bins)) if nch[i] >= max(nch)*min_bw_frac])

        # For well-sampled data (2022 onward) with very low fluence, perform some sigma-clipping to reduce the amount of noise being multiplied in
        if fluence*scale_factor <  0.5 and int(obsid) > 1340639535:
            fl = np.array([fluence_bins[i] for i in range(len(fluence_bins)) if freq_range[i] >= max(freq_range)*min_bw_frac])
            rms = np.nanstd(sigma_clip(fl))
            fluence = np.sum(fl[fl > 3*rms])

        row = {
            'num_obs': 1,
            'pulse_number': new_pulse_number,
            'utcs': [utc],
            'toa': f"{toa.mjd:.5f}",
            'telescopes': [telescope],
            'midfreq': [f"{midfreq:.0f}"],
            'peak': f"{round_to_sf(peak_flux_density, sf=sf)}",
            'peak_1GHz': f"{round_to_sf(peak_flux_density*scale_factor.value*1e3, sf=sf)}", # in mJy
            'fluence': f"{round_to_sf(fluence, sf=sf)}",
            'fluence_1GHz': f"{round_to_sf(fluence*scale_factor.value, sf=sf)}",
        }

        if SHOW_PLOTS:
            plt.plot(t, lc, label=f"LC #1 (dt = {dt})")
            #plt.plot(t, freq_range, label="Freq range #1 (MHz)")
            #plt.legend()
            #plt.show()
    else:
        # Pop the previous row from the table
        row = table.pop()
        fluence_bins_to_plot.pop()

        # Update row values
        row['num_obs'] += 1
        row['utcs'].append(utc)
        row['telescopes'].append(telescope)
        row['midfreq'].append(f"{midfreq:.0f}")

        # For the peak and fluence, we have to update the lightcurve first

        # Combine the light curve with the previous one
        # "Combine" means a weighted sum (equivalent to forming the light curve from both dynamic spectra)
        new_t = lightcurve_data[:,0]
        new_lc = lightcurve_data[:,1]
        new_nch = lightcurve_data[:,2] # The number of channels that contributed to each point in the lightcurve
        new_dt = new_t[1] - new_t[0]
        new_bw = params['bw']
        new_freq_range = new_nch*new_bw
        new_fluence_bins = new_lc * new_dt

        for new_i in range(len(new_t)): # Loop through the new time axis
            t_idx = round((new_t[new_i] - t[0])/dt) # The bin number relative to the previous lightcurve
            if t_idx in t_idxs: # If this bin overlaps with the previous lightcurve
                i = t_idxs.index(t_idx)
                lc[i] = (lc[i]*freq_range[i]*dt + new_lc[new_i]*new_freq_range[new_i]*new_dt) / (freq_range[i]*dt + new_freq_range[new_i]*new_dt) # Weighted sum
                fluence_bins[i] = lc[i] * dt
                freq_range[i] += new_freq_range[new_i] # Update the channel count
            else: # If this bin does not overlap with the previous lightcurve, then none of the others will either, so just append the rest in bulk
                t = np.hstack((t, new_t[new_i:]))
                lc = np.hstack((lc, new_lc[new_i:]))
                freq_range = np.hstack((freq_range, new_freq_range[new_i:]))
                fluence_bins = np.hstack((fluence_bins, new_fluence_bins[new_i:]))
                break

        # Same peak and fluence calculations as before
        peak_flux_density = np.max(lc)
        fluence = np.sum([fluence_bins[i] for i in range(len(fluence_bins)) if freq_range[i] >= max(freq_range)*min_bw_frac])

        row['peak'] = f"{round_to_sf(peak_flux_density, sf=sf)}"
        row['peak_1GHz'] = f"{round_to_sf(peak_flux_density*scale_factor.value*1e3, sf=sf)}" # in mJy
        row['fluence'] = f"{round_to_sf(fluence, sf=sf)}"
        row['fluence_1GHz'] = f"{round_to_sf(fluence*scale_factor.value, sf=sf)}"


        if SHOW_PLOTS:
            plt.plot(new_t, new_lc, label=f"LC #2 (dt = {new_dt})")
            #plt.plot(new_t, new_freq_range, label="Freq range #2 (MHz)")
            #plt.legend()
            #plt.show()

            plt.plot(t, lc, label=f"Combined LC (dt = {dt})")
            #plt.plot(t, freq_range, label="Combined freq range (MHz)")

            plt.xlabel("Time (GPS second)")
            #plt.ylabel("Flux density (Jy) & Freq range (MHz)")
            plt.ylabel("Flux density (Jy)")
            plt.legend()
            plt.show()

    # Add the row to the table
    table.append(row)
    prev_pulse_number = new_pulse_number
    fluence_bins_to_plot.append({
        "obsid": obsid,
        "utc": utc,
        "fluence_bins": fluence_bins,
        "freq_range": freq_range,
    })

# Write out 'fluence lightcurves', just to have something to eyeball to check for any irregularities
print("Plotting fluence bins...")
for fluence_bin_dict in fluence_bins_to_plot:

    # Pull out shorthand variables from dictionary
    obsid = fluence_bin_dict['obsid']
    utc = fluence_bin_dict['utc']
    fluence_bins = fluence_bin_dict['fluence_bins']
    freq_range = fluence_bin_dict['freq_range']

    time_bins = np.arange(len(fluence_bins))
    included = freq_range > max(freq_range)*min_bw_frac

    fig, ax = plt.subplots(1,1)
    plt.plot(time_bins, fluence_bins, label="Fluence in each bin")
    plt.plot(time_bins[included], fluence_bins[included], label="Included in fluence calculation")
    plt.xlabel("Time bin number")
    plt.ylabel("Fluence (Jy s)")
    plt.title(f"{utc}")
    plt.legend()
    plt.savefig(f'{obsid}_fluencebins.png')
    plt.close()

# Format the table to LaTeX table format
print("Writing pulse_table.tex...")
with open("pulse_table.tex", "w") as tex:

    # The column headers:
    print("Pulse & UTC & Barycentered & Telescope & Frequency & \\multicolumn{2}{c}{Peak flux density} & \\multicolumn{2}{c}{Fluence} \\\\", file=tex)

    # The column subheaders
    print("number & & TOA & & & At freq & At 1 GHz & At freq & At 1 GHz \\\\", file=tex)
    print(" & & (MJD) & & (MHz) & (Jy) & (mJy) & (Jy s) & (Jy s) \\\\", file=tex)

    print("\\hline", file=tex)

    # The data
    def multirow(val):
        return f"\\multirow{{{row['num_obs']}}}{{*}}{{{val}}}"

    for row in table:

        if row['num_obs'] == 1:
            print(f"{row['pulse_number']} & {row['utcs'][0]} & {row['toa']} & {row['telescopes'][0]} & {row['midfreq'][0]} & {row['peak']} & {row['peak_1GHz']} & {row['fluence']} & {row['fluence_1GHz']} \\\\", file=tex)
        else:
            print(f"{multirow(row['pulse_number'])} & {row['utcs'][0]} & {multirow(row['toa'])} & {row['telescopes'][0]} & {row['midfreq'][0]} & {multirow(row['peak'])} & {multirow(row['peak_1GHz'])} & {multirow(row['fluence'])} & {multirow(row['fluence_1GHz'])} \\\\", file=tex)
            for i in range(1, row['num_obs']):
                print(f"  & {row['utcs'][i]} & & {row['telescopes'][i]} & {row['midfreq'][i]} & & & & \\\\", file=tex)


# Also, make a simple CSV
print("Writing pulse_table.csv...")
with open("pulse_table.csv", "w") as csv:
    # Write the header line
    print("Pulse number,UTC,Barycentred TOA (MJD),Telescope,Frequency (MHz),Peak flux density at freq (Jy),Peak flux density at 1 GHz (mJy),Fluence at freq (Jy s),Fluence at 1 GHz (Jy s)", file=csv)

    # Write the rows
    for row in table:
        for i in range(row['num_obs']):
            print(f"{row['pulse_number']},{row['utcs'][i]},{row['toa']},{row['telescopes'][i]},{row['midfreq'][i]},{row['peak']},{row['peak_1GHz']},{row['fluence']},{row['fluence_1GHz']}", file=csv)
