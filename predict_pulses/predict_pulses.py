import matplotlib as mpl

mpl.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

import numpy as np
import pandas as pd

# import math

from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation

from glob import glob

import urllib.request
import json

import logging
import argparse

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(module)s:%(levelname)s:%(lineno)d %(message)s")
logger.setLevel(logging.INFO)

BASEURL = "http://ws.mwatelescope.org/"

# MWA location from CONV2UVFITS/convutils.h
LAT = -26.703319
LON = 116.67081
ALT = 377.0

LOCATION = EarthLocation.from_geodetic(
    lat=LAT * u.deg, lon=LON * u.deg, height=ALT * u.m
)


def getmeta(servicetype="metadata", service="obs", params=None):
    """Given a JSON web servicetype ('observation' or 'metadata'), a service name (eg 'obs', find, or 'con')
    and a set of parameters as a Python dictionary, return a Python dictionary containing the result.
    """
    if params:
        # Turn the dictionary into a string with encoded 'name=value' pairs
        data = urllib.parse.urlencode(params)
    else:
        data = ""

    # Get the data
    try:
        result = json.load(
            urllib.request.urlopen(BASEURL + servicetype + "/" + service + "?" + data)
        )
    except urllib.error.HTTPError as err:
        print(
            "HTTP error from server: code=%d, response:\n %s" % (err.code, err.read())
        )
        return
    except urllib.error.URLError as err:
        print("URL or network error: %s" % err.reason)
        return

    # Return the result dictionary
    return result


def barycentre_time_delta(src_coord, location, time):
    temp_time = Time(time, location=location)
    ltt_bary = temp_time.light_travel_time(src_coord)

    return ltt_bary


def predict_times(
    ra=279.75833333,
    dec=-10.53041667,
    separation=50.0,
    starttime=1199000000,
    stoptime=1210000000,
    minchan=60,
    maxchan=180,
    period=21.97,
    pulseref=59780.80241087963,
    pulserefgps=None,
    obslength=596,
    output="periodic_search.txt",
    vcs=False,
):

    # Transient location
    c1 = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="fk5")

    # Transient period
    period = period * 60  # seconds

    # Reference pulse start time
    if pulserefgps is None:
        t0 = Time(pulseref, format="mjd")
    else:
        t0 = Time(pulserefgps, format="gps")

    logger.info(f"Period estimate: {period/60.:3.4f} minutes")
    logger.info(f"Reference pulse time (MJD) {t0.mjd:8.4f}")
    logger.info(f"Reference pulse time (GPS) {t0.gps:8.4f}")

    records = []

    assert starttime < stoptime, "Startime needs to be before stoptime"

    minimum_n = np.floor((starttime - t0.gps) / period).astype(int)
    maximum_n = np.ceil((stoptime - t0.gps) / period).astype(int)

    obslist = []

    logger.info(f"Minimum Pulse Count {minimum_n}")
    logger.info(f"Maximum Pulse Count {maximum_n}")

    # Get the list of obsids that potentially are involved
    for n in range(minimum_n, maximum_n):
        predicted_time = t0.gps + n * period
        predicted_time = Time(predicted_time, format="gps")

        logger.debug(f"{n=}, {predicted_time.gps=}")

        try:
            olist = getmeta(
                service="find",
                params={
                    "mintime": int(predicted_time.gps - obslength),
                    "maxtime": int(predicted_time.gps + obslength),
                    "dict": 1,
                },
            )

            logger.debug(f"Returned {len(olist)=}")

            # Drop out the time here
            obslist.extend(
                [
                    (obs, predicted_time)
                    for obs in olist
                    if obs["mwas.starttime"] < predicted_time.gps < obs["mwas.stoptime"]
                ]
            )
            logger.info(f"Observation list now {len(obslist)=}")
        except:
            pass

    # Make sure we have something to work with
    if obslist is None or len(obslist) == 0:
        print("No results returned to cycle through. ")
        return

    # Now the tests
    for (obs, predicted_time) in obslist:

        # logger.debug(
        #     obs["mwas.starttime"],
        #     obs["sm.ra_pointing"],
        #     obs["sm.dec_pointing"],
        #     obs["rfs.frequencies"][12],
        #     obs["numfiles"],
        #     obs["mwas.mode"],
        # )

        # Do the common checks
        freq_check = minchan <= obs["rfs.frequencies"][12] <= maxchan
        no_f_check = obs["numfiles"] > 2

        logger.debug(f"{obs['mwas.starttime']=} {freq_check=} {no_f_check}=")

        # Actually check the common cases
        if False in (freq_check, no_f_check):
            continue

        if vcs is True:
            # VCS search
            logger.debug(f"VCS checks {obs['mwas.mode']=}")
            if not obs["mwas.mode"] in ("VOLTAGE_START", "MWAX_VCS"):
                logger.debug(f"VCS check failed")
                continue

        else:
            logger.debug(
                f"Continuum mode checks {obs['mwas.mode']=} {obs['mwas.dataquality']=}"
            )
            if not (obs["mwas.mode"] == "HW_LFILES" and obs["mwas.dataquality"] < 4):
                continue

                # Some of the observations are antenna tests or VCS or drift scans so they don't have RAs and Decs
                #                    if obs['sm.ra_pointing'] is not None and obs['sm.dec_pointing'] is not None:
                #                        if obs['mwas.ra_phase_center'] is not None and obs['mwas.dec_phase_center'] is not None:
                #                            c2 = SkyCoord(obs['mwas.ra_phase_center'], obs['mwas.dec_phase_center'], unit = (u.deg, u.deg), frame='fk5')
                # There can be only one observation that meets these criteria so break the loop if found

        c2 = SkyCoord(
            obs["sm.ra_pointing"],
            obs["sm.dec_pointing"],
            unit=(u.deg, u.deg),
            frame="fk5",
        )

        sep = c1.separation(c2)
        logger.debug(f"Sepeation {sep.deg=} degrees")

        if sep.deg >= separation:
            continue

        ltt_bary = barycentre_time_delta(c1, LOCATION, predicted_time)
        logger.debug(f"Refined light-travel-time is {ltt_bary.value=}")

        predicted_time_bary = predicted_time + ltt_bary

        # If we make it to here, we save the record
        records.append(
            dict(
                obsid=obs["mwas.starttime"],
                centchan=obs["rfs.frequencies"][12],
                separation=sep.deg,
                mjd_pulse=predicted_time.mjd,
                gps_pulse=predicted_time.gps,
                bary_mjd_pulse=predicted_time_bary.mjd,
                bary_gps_pulse=predicted_time_bary.gps,
            )
        )

    if len(records) == 0:
        logger.error("No observations found. ")

        import sys

        sys.exit()

    records_df = pd.DataFrame(records)
    records_df.to_csv(output, sep=" ")

    logger.info(f"Have written dataframe: ")
    logger.info(f"  - {len(records_df)=}")
    logger.info(f"  - {records_df.columns=}")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("Time (GPS)")
    ax.set_ylabel("Separation (deg)")
    ax.scatter(records_df["gps_pulse"], records_df["separation"])
    # ax.legend()
    fig.savefig("predictions.png", bbox_inches="tight")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Input/output files")
    group1.add_argument(
        "--ra",
        dest="ra",
        type=float,
        default=279.75833333,
        help="RA / longitude value (deg)",
    )
    group1.add_argument(
        "--dec",
        dest="dec",
        type=float,
        default=-10.53041667,
        help="Dec / latitude value (deg)",
    )
    group1.add_argument(
        "--starttime",
        dest="starttime",
        type=int,
        default=1199000000,
        help="minimum GPS starttime for search (default = 1199000000)",
    )
    group1.add_argument(
        "--stoptime",
        dest="stoptime",
        type=int,
        default=1210000000,
        help="minimum GPS stoptime for search (default = 1210000000)",
    )
    group1.add_argument(
        "--minchan",
        dest="minchan",
        type=int,
        default=60,
        help="minimum coarse channel number (default = 60)",
    )
    group1.add_argument(
        "--maxchan",
        dest="maxchan",
        type=int,
        default=180,
        help="maximum coarse channel for search (default = 180)",
    )
    group1.add_argument(
        "--obslength",
        dest="obslength",
        type=int,
        default=596,
        help="maximum observation length to search (smaller is faster but might miss detections) (default = 596)",
    )
    group1.add_argument(
        "--period",
        dest="period",
        type=float,
        default=21.97,
        help="Period in minutes (default = 18.18480809)",
    )
    group1.add_argument(
        "--pulseref",
        dest="pulseref",
        type=float,
        default=59780.80241087963,
        help="Reference for pulse arrival time in MJD (default = 58151.9593)",
    )
    group1.add_argument(
        "--pulserefgps",
        dest="pulserefgps",
        type=float,
        default=None,
        help="Reference for pulse arrival time in MJD (default = unused, overwrites pulseref if set)",
    )
    group1.add_argument(
        "--separation",
        dest="separation",
        type=float,
        default=50.0,
        help="Maximum allowable sky separation between pointing centre and source (default = 50deg)",
    )
    group1.add_argument(
        "--vcs",
        dest="vcs",
        action="store_true",
        default=False,
        help="Search for VCS observations rather than HW_LFILES",
    )
    group1.add_argument(
        "--output",
        dest="output",
        type=str,
        default="periodic_search.txt",
        help="Output text file for search results (default = periodic_search.txt",
    )
    group1.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help="Enable more logging",
    )

    results = parser.parse_args()

    if results.verbose:
        logger.info("Setting to debug mode")
        logger.setLevel(logging.DEBUG)

    predict_times(
        ra=results.ra,
        dec=results.dec,
        separation=results.separation,
        starttime=results.starttime,
        stoptime=results.stoptime,
        minchan=results.minchan,
        maxchan=results.maxchan,
        period=results.period,
        pulseref=results.pulseref,
        pulserefgps=results.pulserefgps,
        obslength=results.obslength,
        output=results.output,
        vcs=results.vcs,
    )
