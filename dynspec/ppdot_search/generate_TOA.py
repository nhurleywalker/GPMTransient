import numpy as np
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

    # Generate a Gaussian with the specified width, but with the
    #width_

if __name__ == '__main__':
    # Parse the command line
    parser = argparse.ArgumentParser(description='Generate TOAs from lightcurves using a basic low-pass Gaussian filter')
    parser.add_argument('yaml', type=argparse.FileType('r'), help='The yaml file containing the meta information for the observation')
    parser.add_argument('lightcurve', type=str, help='The two-column ASCII file containing the lightcurve data')
    parser.add_argument('kernel_width_in_sec', type=float, help='The 1σ width of the Gaussian used for the filter (in seconds)')

    args = parser.parse_args()

    main(args)
