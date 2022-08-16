#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from scipy.interpolate import interp1d
from scipy.optimize import fmin
from glob import glob
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
import multiprocessing
import copy

DMP = True
try:
    import pygedm
except ImportError:
    DMP = False
    print("No position-dependent DM trials available")

# The DM delay formula assumes time in seconds, frequencies in MHz
def calc_dmdelay(DM, flo, fhi):
    return 4.148808e3*DM*(1/(flo*flo) - 1/(fhi*fhi))

def sc(data):
    '''
    Perform some sigma-clipping on a light curve
    '''
    std = np.std(data)
    std = np.std(data[np.abs(data) < 3*std])
    return std

class Dynspec:
    def __init__(self, dynspec, sample_time, freqlo, bw, time_offset=0, freq_offset=0, dm=0):
        self.TIMEAXIS=1
        self.FREQAXIS=0
        self.dynspec=dynspec.T
        self.Nf = self.dynspec.shape[self.FREQAXIS]
        self.Nt = self.dynspec.shape[self.TIMEAXIS]
        self.create_time_axis(sample_time, offset=time_offset)
        self.create_freq_axis(bw, freqlo=freqlo)
        self.dm = dm # The initial DM of the loaded data
        self.fscrunched = None

    def create_time_axis(self, sample_time, offset=0):
        self.t = np.arange(self.Nt)*sample_time + offset + sample_time/2
        self.dt = sample_time

    def create_freq_axis(self, bw, freqlo=0):
        self.f = np.arange(self.Nf)*bw + freqlo
        self.df = bw

    def plot(self, ax):
        extent = (
                self.t[0]  - self.dt/2,
                self.t[-1] + self.dt/2,
                self.f[0]  - self.df/2,
                self.f[-1] + self.df/2)
        ax.imshow(self.dynspec, aspect='auto', interpolation='none', origin='lower', extent=extent)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Frequency (MHz)')

    def plot_lightcurve(self, ax):
        if self.fscrunched is not None:
            ax.plot(self.t, self.fscrunched)

    def set_freq_ref(self, freq_ref):
        if freq_ref is None or freq_ref == 'centre':
            self.freq_ref = np.mean(self.f)
        elif freq_ref == 'low':
            self.freq_ref = self.f[0]
        elif freq_ref == 'high':
            self.freq_ref = self.f[-1]
        else:
            try:
                self.freq_ref = float(freq_ref)
            except:
                meanfreq = np.mean(self.f)
                print("Warning: could not interpret {} as frequency, setting the reference frequency to {} MHz".format(freq_ref, meanfreq))
                self.freq_ref = meanfreq

    def calc_dm_shifts(self, dm, freq_ref=None):
        # Set the reference frequency
        self.set_freq_ref(freq_ref)

        # Calculate how much the DM is to change by
        dm_diff = dm - self.dm

        # The DM delay formula assumes time in seconds, frequencies in MHz
        return calc_dmdelay(dm_diff, self.freq_ref, self.f)

    def dedisperse(self, dm, freq_ref=None):
        shifts = self.calc_dm_shifts(dm, freq_ref=freq_ref)

        # FFT along the time axis
        ffted = np.fft.rfft(self.dynspec, axis=self.TIMEAXIS)

        # Calculate phase ramp equivalent to a time shift
        fftfreq = np.fft.rfftfreq(self.Nt, d=self.dt)
        F, T = np.meshgrid(fftfreq, shifts)
        phase_ramp = np.exp(-2j*np.pi*F*T)
        #print(T[0,:], F[0,:], phase_ramp[0,:])

        # Apply the phase ramp
        ffted *= phase_ramp
        #print(ffted[0,:])

        # Inverse FFT to get the dedispersed dynamic spectrum
        self.dynspec = np.fft.irfft(ffted, n=self.Nt, axis=self.TIMEAXIS)

        # Record the current DM
        self.dm = dm

    def fscrunch(self):
        self.fscrunched = np.mean(self.dynspec, axis=self.FREQAXIS)

    def get_time_at_infinite_frequency(self):
        self.dmdelay = calc_dmdelay(self.dm, self.freq_ref, np.inf)
        return self.t - self.dmdelay

class DMCurve():
    def __init__(self, dynspec):
        self.dynspec = dynspec

    def run_dmtrials(self, dmtrials, freq_ref=None):
        self.dms = np.arange(*dmtrials)
        self.peak_snrs = []
        for dm in self.dms:
            self.dynspec.dedisperse(dm, freq_ref=freq_ref)
            self.dynspec.fscrunch()
            self.peak_snrs.append(np.max(self.dynspec.fscrunched))
        self.peak_snrs = np.array(self.peak_snrs)

    def calc_best_dm(self):
        '''
        Assumes run_dmtrials has already been called
        '''
        # Interpolate the DM curve and get the max value
#        try:
#        interped = interp1d(self.dms, self.peak_snrs, kind='cubic')
#        closest_best_dm_idx = np.argmax(self.peak_snrs, keepdims=True)[0]
# Old numpy compatibility
        closest_best_dm_idx = np.argmax(self.peak_snrs)
#        self.best_dm = fmin(lambda dm: -interped(dm), self.dms[closest_best_dm_idx])
        self.best_dm = self.dms[closest_best_dm_idx]
        self.best_snr = self.peak_snrs[closest_best_dm_idx]
#        except ValueError:
#        self.best_dm = [0.0]

def _main(pargs):
    """
    A shallow wrapper for main

    Parameters
    ----------
    args : list
        A list of arguments for main

    Returns
    -------
    None
    """
    # an easier to debug traceback when multiprocessing
    # thanks to https://stackoverflow.com/a/16618842/1710603
    try:
        return main(*pargs)
    except:
        import traceback
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

def main(n, coords, dynspec, args):
    # Load the data into a "dynspec" object
    dynspec = Dynspec(dynspec,
            sample_time=args.sample_time,
            freqlo=args.freqlo,
            bw=args.bw,
            time_offset=args.t0)

    # Run a bunch of DM trials to create a DM curve
    dmcurve = DMCurve(dynspec)
    print(args.dms)
    dmcurve.run_dmtrials(args.dms, freq_ref=args.freq_ref)
    dmcurve.calc_best_dm()
    DM = dmcurve.best_dm
    peak = dmcurve.best_snr

    # Dedisperse the spectrum
    dynspec.dedisperse(DM, freq_ref=args.freq_ref)
    dynspec.fscrunch()

    snr = peak / sc(dynspec.fscrunched)

    if args.no_plots == False and snr > 6:
        # Plot the DM curve
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.plot(dmcurve.dms, dmcurve.peak_snrs)
        ax.set_xlabel("DM (pc/cm^3)")
        ax.set_ylabel("Peak flux density (a.u.)")
        fig.savefig(f"lc_best_DM_{n:04d}.png", bbox_inches="tight")

        # Plot the dynamic spectrum at the given/best DM
        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
        dynspec.plot_lightcurve(axs[0])
        dynspec.plot(axs[1])
        fig.suptitle(f'x, y = {coords[1], coords[0]}, DM = {DM:.1f} pc/cm^3')
        axs[0].set_yticks([])
        fig.savefig(f"ds_best_DM_{n:04d}.png", bbox_inches="tight")

    # If requested, write out a time series of the frequency-scrunched
    # lightcurve
    if args.lightcurve is not None and snr > 6:
        # Get the time of the first bin referenced to infinite frequency
        timeaxis = dynspec.get_time_at_infinite_frequency()
        timeaxis += args.bc_corr # Add the barycentric correction

        # Create verbose header for lightcurve output files
        header = "Created with:\n"
        header += '  ' + ' '.join(sys.argv) + '\n\n'
        header += 'This time series has been dedispersed to {} pc/cm^3\n'.format(DM)
        header += 'Using barycentric correction of {} s\n'.format(args.bc_corr)
        header += 'Using dedispersion delay of {} s (for reference frequency {})\n\n'.format(dynspec.dmdelay, dynspec.freq_ref)
        header += "Time (s) | Flux density (a.u.)"

        # Construct array to be written out and write it out
        lightcurve = np.array([timeaxis, dynspec.fscrunched]).T
        np.savetxt(args.lightcurve, lightcurve, header=header)

    if args.output is not None:
        np.savetxt(args.output, dynspec.dynspec)

    return(n, coords, DM, peak)

if __name__ == "__main__":
    # Parse the command line
    parser = argparse.ArgumentParser(description='Dedisperse a set of images where each pixel is a dynamic spectrum')
    parser.add_argument('--dms', type=float, nargs='*', help='DM trials (pc/cm^3) [[start=0], stop, [step=1]] (i.e. same argument structure as s NumPy''s arange function). If a single value is given, the dynamic spectrum is dedispersed to this DM and no search is performed')
    parser.add_argument('--sample_time', type=float, default=0.5, help='The time of one sample (s)')
    parser.add_argument('--freqlo', type=float, default=139.52, help='The centre frequency of the lowest channel (MHz)')
    parser.add_argument('--freq_ref', type=str, default='centre', help='The reference frequency used during dedispersion. Either a frequency in MHz, or one of [\'low\', \'centre\', \'high\'] (default=\'centre\')')
    parser.add_argument('--bw', type=float, default=1.28, help='The channel width (MHz)')
    parser.add_argument('--output', type=argparse.FileType('w'), help='The file to which the dedispersed dynamic spectrum will be written')
    parser.add_argument('--input', type=str, help='The file stem to glob on; expecting a format like <stem>-t????-????-image.fits (i.e. WSClean output format)')
    parser.add_argument('--lightcurve', type=argparse.FileType('w'), help='Write out the frequency-scrunched, dispersion-corrected lightcurve to the named file. "Dispersion-corrected" means using infinite frequency as reference')
    parser.add_argument('--t0', type=float, default=0, help='The left (early) edge of the first time bin')
    parser.add_argument('--no_plots', action='store_true', help='Do NOT make Matplotlib plots of the DM curve and dedispersed spectrum')
    parser.add_argument('--bc_corr', type=float, default=0, help='Barycentric correction to apply (in seconds)')
    parser.add_argument('--cores', dest='cores', default=None, type=int,
                        help="Number of cores to use (default = autodetect")

    args = parser.parse_args()
    if args.cores is None:
        cores = multiprocessing.cpu_count()
    else:
        cores = args.cores
    print(f"Using {cores} cores...")

# Make a copy of the default DMs
    dms = args.dms

    # read in the images
    hdus = sorted(glob(f"{args.input}-t????-????-image.fits"))
    h = hdus[-1]
    st = h.split("-t")[1]
    tmax = int(st.split("-")[0])
    cmax = int(st.split("-")[1])
    a = fits.open(h)
    xmax = a[0].data.shape[2]
    ymax = a[0].data.shape[3]
    w = WCS(a[0].header, naxis=2)

    # make a (hyper-) cube
    arr = np.empty((tmax, cmax, ymax, xmax))
    print(arr.shape)
    #for t in range(0, tmax):
    #    for c in range(0, cmax):
# Fast testing
    for t in range(0, 3):
        for c in range(0, 3):
            h = fits.open(f"{args.input}-t{t:04d}-{c:04d}-image.fits")
            arr[t, c, :, :] = h[0].data

    # When using a pool, you need to send the arguments as a tuple through the
    # wrapping function, so we set up the arguments here before running the pool
    pargs = []
    # We need to order our arguments by an index, since the results will be returned
    # in whatever order the pooled tasks finish
    n = 0
#    for x in range(0,xmax):
#        for y in range(0,ymax):
# Fast testing
    for x in range(499,500):
        for y in range(0,1000):
            coords = w.pixel_to_world(y, x)
            if np.abs(coords.galactic.b.value) < 2.1:
                if DMP is True:
                    DM, tau_sc = pygedm.dist_to_dm(coords.galactic.l.value, coords.galactic.b.value, 10000, method='ymw16')
                    xargs = copy.deepcopy(args)
                    xargs.dms[1] = DM.value
                pargs.append((n, [x,y], arr[:, :, x, y], xargs))
                n+=1

    # make a multiprocessing pool
    # start a new process for each task
    pool = multiprocessing.Pool(processes=cores, maxtasksperchild=1)
    try:
        # chunksize=1 ensures that we only send a single task to each process
        results = pool.map_async(_main, pargs, chunksize=1).get(timeout=10000000)
    except KeyboardInterrupt:
        pool.close()
        sys.exit(1)
    pool.close()
    pool.join()

    indices, coords, DMs, peaks = map(list, zip(*results))
    # Order correctly
    ind = np.argsort(indices)
    DMs = np.array(DMs)
    DMs = DMs[ind]
    peaks = np.array(peaks)
    peaks = peaks[ind]

    a[0].data = np.zeros(a[0].data.shape)

    for c, DM in zip(coords, DMs):
        a[0].data[:,:,c[0],c[1]] = DM
    a.writeto("DM_map.fits", overwrite=True)

    for c, p in zip(coords, peaks):
        a[0].data[:,:,c[0],c[1]] = p
    a.writeto("SNR_map.fits", overwrite=True)

