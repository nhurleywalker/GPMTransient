import subprocess
import argparse
from astropy.time import Time

def gps2mjd(gpstimes):
    return [Time(gpstime, format='gps').mjd for gpstime in gpstimes]

def bc_corr(ra_hrs, dec_degs, mjds, ephemeris_file, path_to_binary='bc_corr'):
    '''
    This function uses an external program which lives in
    https://github.com/robotopia/natashas_mystery_source
    '''
    return [float(subprocess.run([path_to_binary,
        str(ra_hrs),
        str(dec_degs),
        str(mjd),
        ephemeris_file], capture_output=True).stdout) for mjd in mjds]

if __name__ == "__main__":
    # Parse the command line
    parser = argparse.ArgumentParser(description='Calculate the barycentric correction towards a given source at specific times')
    parser.add_argument('ra', type=float, help='The RA in decimal hours')
    parser.add_argument('dec', type=float, help='The declination in decimal degrees')
    parser.add_argument('ephemeris', type=str, help='The ephemeris file (e.g. http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp)')
    parser.add_argument('gpstimes', type=float, nargs='*', help='The GPS times at which the barycentric correction is to be calculated')
    parser.add_argument('--path_to_bc_corr', type=str, default='bc_corr', help='The path to the bc_corr binary which does the calculation')

    args = parser.parse_args()

    mjds = gps2mjd(args.gpstimes)

    corrected = bc_corr(args.ra, args.dec, mjds, args.ephemeris, path_to_binary=args.path_to_bc_corr)

    for c in corrected:
        print(c)
