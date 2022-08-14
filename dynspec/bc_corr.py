import numpy as np
import argparse
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.constants import c
from spiceypy.spiceypy import spkezr, furnsh, j2000, spd, unload

def bc_corr(coord, times, ephemeris_file='de430.bsp'):
    '''
    coord - SkyCoord object (from astropy) representing the location of the source
    times - Time object (from astropy)
    ephemeris_file - 
    '''
    try:
        furnsh(ephemeris_file)
    except:
        raise Exception("Cannot load ephemeris file {}\n".format(ephemeris_file))
    jds = np.array([time.mjd + 2400000.5 for time in times])
    ets = (jds - j2000())*spd()
    r_earth = [spkezr("earth", et, "j2000", "NONE", "solar system barycenter")[0][:3] for et in ets]
    r_src_normalised = [np.cos(coord.ra.rad)*np.cos(coord.dec.rad),
            np.sin(coord.ra.rad)*np.cos(coord.dec.rad),
            np.sin(coord.dec.rad)]
    delays = [np.dot(r_earth, r_src_normalised) * 1e3 * u.meter / c] # (1e3 = km->m)
    return delays

if __name__ == "__main__":
    # Parse the command line
    parser = argparse.ArgumentParser(description='Calculate the barycentric correction towards a given source at specific times')
    parser.add_argument('ra', type=float, help='The RA in decimal hours')
    parser.add_argument('dec', type=float, help='The declination in decimal degrees')
    parser.add_argument('--ephemeris', type=str, default='de430.bsp', help='The ephemeris file (e.g. http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp)')
    parser.add_argument('gpstimes', type=float, nargs='*', help='The GPS times at which the barycentric correction is to be calculated')

    args = parser.parse_args()

    times = [Time(gpstime, format='gps') for gpstime in args.gpstimes]

    coord = SkyCoord(ra=args.ra*u.hr, dec=args.dec*u.deg, frame='icrs')
    corrected = bc_corr(coord, times, ephemeris_file=args.ephemeris)

    print(corrected)
