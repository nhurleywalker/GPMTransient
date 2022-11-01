import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import yaml
import argparse
from astropy.time import Time

import sys
from os import path
sys.path.append("..")
import dedisperse_dynspec as dd

tempo_telescope_names = {"murriyang": "pks",
                         "parkes": "pks",
                         "pks": "pks",
                         "mwa": "mwa",
                         "meerkat": "meerkat",
                         "vla": "vla",
                         "gmrt": "gmrt"}

def main(args):
    # Read in the metadata from the YAML file
    metadata = dd.parse_yaml(args.yaml)

    # Generate a Gaussian with the specified width, with the same time resolution as the data
    dt = metadata['sample_time']
    nsamples_gaussian_half = int(np.ceil(5*args.kernel_width/dt))

    kernel_t = np.arange(-nsamples_gaussian_half, nsamples_gaussian_half+1)*dt
    kernel_t_rolled = np.roll(kernel_t, -nsamples_gaussian_half)
    kernel = norm.pdf(kernel_t_rolled, scale=args.kernel_width)

    # Read in the lightcurve and FFT it
    lightcurve_data = np.loadtxt(args.lightcurve)
    lightcurve = lightcurve_data[:,1]
    lightcurve_rffted = np.fft.rfft(lightcurve, norm='ortho')

    # Convert the kernel into a filter
    kernel_rffted = np.fft.rfft(kernel, norm='ortho')

    # If the lightcurve has more samples than the kernel, zero pad its FFT. Otherwise, truncate it.
    Nk = len(kernel_rffted)
    Nl = len(lightcurve_rffted)
    if Nl > Nk:
        lightcurve_rffted = lightcurve_rffted[:Nk]
    elif Nk > Nl:
        kernel_rffted = kernel_rffted[:Nl]

    # Apply the filter
    lightcurve_filtered_rffted = lightcurve_rffted * kernel_rffted

    # If the precision wanted is higher than the data time resolution,
    # calculate the number of time bins in the result
    N = len(lightcurve)
    newN = N
    if args.toa_precision is not None:
        if args.toa_precision < dt:
            newN *= int(np.round(dt / args.toa_precision))
    new_dt = dt*(N/newN)

    # Inverse FFT
    lightcurve_filtered = np.fft.irfft(lightcurve_filtered_rffted, n=newN, norm='ortho')

    # Construct the relevant time axes
    t = lightcurve_data[:,0]
    new_t0 = t[0] - dt/2 + new_dt/2
    new_t = np.arange(newN)*new_dt + new_t0

    # Get the maximum in the filtered lightcurve and call that the TOA
    TOA_idx = np.argmax(lightcurve_filtered)
    TOA_gps = new_t[TOA_idx]

    # Convert the (GPS) TOA to MJD
    TOA_mjd = Time(TOA_gps, format='gps').mjd

    # Convert the error (the width of the kernel) to ms
    TOA_err_us = args.kernel_width*1e6

    # Get the reference frequency (cf. Dynspec.set_freq_ref())
    dynspec = dd.Dynspec(**metadata)
    dynspec.set_freq_ref(metadata['freq_ref'])
    freq_MHz = dynspec.freq_ref

    # Print out the TOA in the format expected of timing software
    print(args.lightcurve, freq_MHz, TOA_mjd, TOA_err_us, tempo_telescope_names[metadata['telescope'].lower()])

    # Make a plot, if requested
    if args.save_plot:
        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
        axs[0].plot(t, lightcurve, label="Original")
        axs[1].plot(new_t, lightcurve_filtered, label="Filtered")
        axs[1].scatter(TOA_gps, lightcurve_filtered[TOA_idx], color='k', label="TOA")
        axs[1].set_xlabel("Time (GPS seconds)")
        axs[0].legend()
        axs[1].legend()

        if args.save_plot == 'SHOW':
            plt.show()
        else:
            plt.savefig(args.save_plot)

if __name__ == '__main__':
    # Parse the command line
    parser = argparse.ArgumentParser(description='Generate TOAs from lightcurves using a basic low-pass Gaussian filter')
    parser.add_argument('yaml', type=argparse.FileType('r'), help='The yaml file containing the meta information for the observation')
    parser.add_argument('lightcurve', type=str, help='The two-column ASCII file containing the lightcurve data')
    parser.add_argument('kernel_width', type=float, help='The 1σ width of the Gaussian used for the filter (in seconds)')
    parser.add_argument('--toa_precision', type=float, default=None, help='The precision of the reported TOA (in seconds). Only changes precision if the desired precision is higher than the data time resolution. Default: Same as data time resolution')
    parser.add_argument('--save_plot', type=str, help='Saves a plot of the original light curve, the filtered light curve, and the TOA. If the argument equals "SHOW", it will only display the plot on the screen.')

    args = parser.parse_args()

    main(args)
