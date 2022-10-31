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
    print(metadata)

    # Generate a Gaussian with the specified width, with the same time resolution as the data
    dt = metadata['sample_time']
    nsamples_gaussian_half = int(np.ceil(5*args.kernel_width_in_sec/dt))

    t = np.arange(-nsamples_gaussian_half, nsamples_gaussian_half+1)*dt
    kernel = norm.pdf(t, scale=args.kernel_width_in_sec)


    # Read in the lightcurve and FFT it
    lightcurve_data = np.loadtxt(args.lightcurve)
    lightcurve = lightcurve_data[:,1]
    lightcurve_rffted = np.fft.rfft(lightcurve)

    # Convert the kernel into a filter
    kernel_rffted = np.fft.rfft(kernel)

    # We want to end up with "dt" time resolution.
    # Therefore, if the lightcurve has more samples than the kernel, zero pad its FFT. Otherwise, truncate it.
    Nk = len(kernel_rffted)
    Nl = len(lightcurve_rffted)
    if Nl > Nk:
        lightcurve_rffted = lightcurve_rffted[:Nk]
    elif Nk > Nl:
        kernel_rffted = kernel_rffted[:Nl]

    # Apply the filter
    lightcurve_filtered = lightcurve_rffted * kernel_rffted
    plt.plot(np.abs(lightcurve_filtered))
    plt.show()

if __name__ == '__main__':
    # Parse the command line
    parser = argparse.ArgumentParser(description='Generate TOAs from lightcurves using a basic low-pass Gaussian filter')
    parser.add_argument('yaml', type=argparse.FileType('r'), help='The yaml file containing the meta information for the observation')
    parser.add_argument('lightcurve', type=str, help='The two-column ASCII file containing the lightcurve data')
    parser.add_argument('kernel_width_in_sec', type=float, help='The 1σ width of the Gaussian used for the filter (in seconds)')

    args = parser.parse_args()

    main(args)
