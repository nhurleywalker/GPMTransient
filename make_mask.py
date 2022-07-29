#!/usr/bin/env python

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy import units as u

import sys

# command line option is a template file, e.g. what you use to make the beam files
inp = sys.argv[1]
mask_file = "{0}_mask.fits".format(inp[0:10])
hdu = fits.open(sys.argv[1])
w = wcs.WCS(hdu[0].header, naxis=2)

# and coordinates of the source-of-interest
coords = sys.argv[2]
loc = w.world_to_pixel(SkyCoord(coords, unit=(u.hour, u.deg), frame="fk5"))

hdu[0].data = np.ones(hdu[0].data.shape, dtype=np.float32)
hdu[0].data[:,:,int(loc[1])-20:int(loc[1])+20,int(loc[0])-20:int(loc[0])+20] = 0.0

hdu.writeto(mask_file, overwrite=True)
