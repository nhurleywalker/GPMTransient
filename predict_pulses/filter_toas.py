import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time, TimeDelta
import argparse

from pint.toa import get_TOAs

# Defined in
# https://github.com/MWATelescope/mwalib/blob/2155c19ad2eebbc6ab2131b012d12e8638d6f493/src/lib.rs#L30
MWA_long = 2.0362898668561042*u.rad
MWA_lat = 0.4660608448386394*u.rad

# Given in top-level README.md
GPM_RA = 18.6505*(15*u.degree)
GPM_Dec = -10.5304*u.degree

def calc_dmdelay(DM, flo, fhi):
    return 4.148808e3*DM*(1/(flo*flo) - 1/(fhi*fhi))

def main():
    # Parse the command line
    parser = argparse.ArgumentParser(description='Filter TOAs according to some criteria.')
    parser.add_argument('timfile', type=argparse.FileType('r'), help='The TOA file (e.g. .tim) containing TOAs to be filtered.')
    parser.add_argument('--HA_range', type=float, nargs=2, help='Keep only TOAs within this hour angle range (MIN, MAX), where MIN and MAX are given in hours, and are in the range [0,24]. To get a range that straddles 0, you can use, e.g., MIN=23, MAX=1.')
    parser.add_argument('--format', choices=['MJD', 'GPS'], default='GPS', help="The output format of the filtered TOAs. Currently, only the values themselves are written out (not the whole TOA line), with support only for MJD or GPS (default) times.")
    parser.add_argument('outfile', help="The file name to write the result to.")

    args = parser.parse_args()

    # Load TOAs
    toas = get_TOAs(args.timfile)
    toa_list = toas.to_TOA_list()
    mjds = np.array([toa.mjd for toa in toa_list])

    if args.HA_range:

        # Setup the MWA's location
        MWA_location = EarthLocation(lat=MWA_lat, lon=MWA_long)

        # Setup the GPM 1839-10's coords
        source_coord = SkyCoord(ra=GPM_RA, dec=GPM_Dec)

        # Convert the source to Alt-Az
        source_altaz = source_coord.transform_to(AltAz(obstime=mjds, location=MWA_location))

        # Get the LHA of the source
        source_LHA = source_altaz.transform_to('fk5').ra.hour - MWA_long.to('degree').value*15

        # Actually apply the hour angle filter
        LHA_min, LHA_max = args.HA_range
        if LHA_min < LHA_max:
            time_mask = np.logical_and(source_LHA >= LHA_min, source_LHA <= LHA_max)
        else:
            time_mask = np.logical_or(source_LHA >= LHA_min, source_LHA <= LHA_max)

        mjds = mjds[time_mask]

    # Print out the result
    with open(args.outfile, "w") as f:
        if args.format == 'GPS':
            for t in mjds:
                print(int(round(t.gps)), file=f)
        elif args.format == 'MJD':
            for t in mjds:
                print(t.mjd, file=f)
        #else: # There shouldn't be an else, because of "choices" in argparse
        #    pass

if __name__ == "__main__":
    main()

