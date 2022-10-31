import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import yaml
import argparse

import sys
from os import path
sys.path.append("..")
import dedisperse_dynspec as dd

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
            newN = N*int(np.round(dt / args.toa_precision))
    new_dt = dt*(N/newN)

    # Inverse FFT
    lightcurve_filtered = np.fft.irfft(lightcurve_filtered_rffted, n=newN, norm='ortho')

    t = lightcurve_data[:,0]
    new_t0 = t[0] - dt/2 + new_dt/2
    new_t = np.arange(new_t0, new_t0 + newN*new_dt, new_dt)
    plt.plot(t, lightcurve, label="Original")
    plt.plot(new_t, lightcurve_filtered*500, label="Filtered")
    plt.legend()
    plt.show()

if __name__ == '__main__':
    # Parse the command line
    parser = argparse.ArgumentParser(description='Generate TOAs from lightcurves using a basic low-pass Gaussian filter')
    parser.add_argument('yaml', type=argparse.FileType('r'), help='The yaml file containing the meta information for the observation')
    parser.add_argument('lightcurve', type=str, help='The two-column ASCII file containing the lightcurve data')
    parser.add_argument('kernel_width', type=float, help='The 1σ width of the Gaussian used for the filter (in seconds)')
    parser.add_argument('--toa_precision', type=float, default=None, help='The precision of the reported TOA (in seconds). Only changes precision if the desired precision is higher than the data time resolution. Default: Same as data time resolution')

    args = parser.parse_args()

    main(args)
