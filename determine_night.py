#!/usr/bin/env python

# Determine the observing night for the follow-up of the Galactic transient

from astropy.time import Time
import argparse

def epoch(obsid=None):
# TODO make the output YYMMDD 
    if obsid is None:
        t = Time.now()
    else:
        t = Time(obsid, format="gps")
# Y2k is 2032
    return(f"{t.ymdhms[0]:04d}-{t.ymdhms[1]:02d}-{t.ymdhms[2]:02d}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--obsid", dest='obsid', type=str, default=None,
                        help="An obsid to determine the epoch of. If not specified, determines epoch of now.")
    args = parser.parse_args()
    print(epoch(args.obsid))
