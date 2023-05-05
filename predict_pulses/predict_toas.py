"""
Predict times of arrival (TOAs) from a given timing model (.par) and a set of known TOAs (.tim).

Strictly speaking, a set of known TOAs are not needed to make TOA predictions because there is
in principle enough information in the model itself. However, this program requires a set of
known TOAs, which serves two purposes:

    1) If the --refit option is given, the known TOAs are used to improve the model before using
       it to predict other TOAs.
    2) For some reason I can't figure out, the phases of the predicted TOAs are (generally) offset
       from the known TOAs by some random amount (constant offset for all predicted TOAs). Without
       a set of known TOAs to compare against, you'd never know there was an offset, so the
       recommended (and unfortunately hacky) solution here is to make use of the --manual_offset
       option to correct for it. The (clunky) workflow is as follows:
          1. generate predicted TOAs without an offset (i.e. don't use --manual_offset)
          2. combine the known and predicted TOAs into a single .tim file
          3. Use TEMPO2/PINT (or whatever) to plot the residuals and use the plots to determine
             what the offset should be.
          4. Re-run this script with --manual_offset.
"""

import numpy as np
import astropy.units as u
from astropy.time import Time, TimeDelta
import argparse

import pint.logging
from pint.fitter import DownhillWLSFitter
from pint.models.model_builder import get_model, get_model_and_toas
from pint.toa import get_TOAs
from pint.simulation import make_fake_toas_fromMJDs
#from pint.solar_system_ephemerides import load_kernel  # <--- I hoped this would speed things up; it didn't

def main():

    # Parse the command line
    parser = argparse.ArgumentParser(description='Predict pulse times within a given range of MJDs for a given ephmeris')
    parser.add_argument('parfile', type=argparse.FileType('r'), help='The ephemeris file (e.g. .par) containing the timing model to use for the prediction.')
    parser.add_argument('timfile', type=argparse.FileType('r'), help='The TOA file (e.g. .tim) containing real TOAs to use as a reference.')
    parser.add_argument('MJDstart', type=float, help='The starting MJD.')
    parser.add_argument('MJDend', type=float, help='The ending MJD.')
    parser.add_argument('--refit', action="store_true", help="Refit the model on the given TOAs before predicting new TOAs.")
    parser.add_argument('--outtimfile', default="out.tim", help="Writing the resulting TOAs to this file (default = 'out.tim')")
    parser.add_argument('--outparfile', help="Writing the (refitted, if --refit is given) timing model to this file. If not supplied, no file is written out.")
    parser.add_argument('--manual_offset', type=float, default=0.0, help="For some unknown reason, PINT generates simulated TOAs that do not necessarily match phase with the given TOAs. Use this flag to add an arbitrary offset (in seconds) to all the simulated TOAs.")
    #parser.add_argument('--planetary_ephemeris', nargs=2, help="The [NAME, PATH] of the planetary ephemeris to use. If not provided, PINT will use its own means for obtaining one -- possibly via very slow network link.")

    args = parser.parse_args()

    '''
    try:
        ephem, path = args.planetary_ephemeris
        load_kernel(ephem, path=path)
    except:
        pass
    '''

    m, t = get_model_and_toas(args.parfile, args.timfile, ephem="DE436")
    #m = get_model(args.parfile)
    #t = get_TOAs(args.timfile)

    if args.refit:
        f = DownhillWLSFitter(t, m)
        f.fit_toas()
        m = f.model

    # Now try to make matching TOAs
    # Initial guess is just pulses with model's F0, and F1=0.
    # These will be the input to PINT's make_fake_toas_fromMJDs() function,
    # which AFAIK "rounds" the given MJDs to the nearest "correct" TOA
    # such to make the residuals zero for the given model
    initial_MJDs = Time(np.arange(args.MJDstart, args.MJDend, (1/m.F0.quantity).to('day').value), format='mjd')
    toas = make_fake_toas_fromMJDs(initial_MJDs, m, freq=150*u.MHz, obs='MWA', include_bipm=True, include_gps=True)
    Δt = TimeDelta(np.full(toas.table["mjd"].shape, args.manual_offset), format='sec')
    toas.adjust_TOAs(Δt)
    toas.write_TOA_file(args.outtimfile)

    if args.outparfile:
        m.write_parfile(args.outparfile)

if __name__ == "__main__":
    main()
