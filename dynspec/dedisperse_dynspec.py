import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from scipy.interpolate import interp1d
from scipy.optimize import fmin, fminbound
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
import yaml
import copy

EPHEMERIS = 'de430.bsp'

# The DM delay formula assumes time in seconds, frequencies in MHz
def calc_dmdelay(DM, flo, fhi):
    return 4.148808e3*DM*(1/(flo*flo) - 1/(fhi*fhi))

class Dynspec:
    def __init__(self, dm=0, **kwargs):
        self.TIMEAXIS=1
        self.FREQAXIS=0

        try:
            transpose = kwargs['transpose']
        except:
            transpose = None

        try:
            dynspec = np.loadtxt(kwargs['input'], ndmin=2)
        except:
            try:
                dynspec = kwargs['dynspec']
            except:
                raise Exception("'input' (filename) or 'dynspec' (2D Numpy array) required for Dynspec initialiser")

        try:
            sample_time = kwargs['sample_time']
        except:
            sample_time = 1.0

        try:
            t0 = kwargs['t0']
        except:
            t0 = 0.0

        try:
            bw = kwargs['bw']
        except:
            bw = 1.0

        try:
            freqlo = kwargs['freqlo']
        except:
            freqlo = 0.0

        self.load_dynspec(dynspec, transpose=transpose)
        self.create_time_axis(sample_time, t0)
        self.create_freq_axis(bw, freqlo)
        self.dm = dm # The initial DM of the loaded data
        self.fscrunched = None

        # Apply any RFI masking
        if 'mask_value' in kwargs:
            mask_value = kwargs['mask_value']
        else:
            mask_value = 0.0

        if 'mask_time_bins' in kwargs:
            self.mask_time_bins(kwargs['mask_time_bins'], mask_value=mask_value)
        if 'mask_freq_bins'in kwargs:
            self.mask_freq_bins(kwargs['mask_freq_bins'], mask_value=mask_value)

        self.freq_ref = None

    def load_dynspec(self, dynspec, transpose=None):
        if transpose is None:
            transpose = False
        # Load the data with NumPy
        self.dynspec = dynspec
        self.mask = np.full(self.dynspec.shape, False, dtype=bool)

        # Swap rows/columns, if specified
        if transpose:
            self.transpose()

        # Record the dimensions of the array
        self.Nf = self.dynspec.shape[self.FREQAXIS]
        self.Nt = self.dynspec.shape[self.TIMEAXIS]

    def mask_time_bins(self, time_bins, mask_value=0.0):
        if time_bins is None:
            return

        if self.TIMEAXIS == 0:
            for idx in time_bins:
                self.dynspec[idx,:] = mask_value
                self.mask[idx,:] = True
        else: # TIMEAXIS == 1
            for idx in time_bins:
                self.dynspec[:,idx] = mask_value
                self.mask[:,idx] = True

    def mask_freq_bins(self, freq_bins, mask_value=0.0):
        if freq_bins is None:
            return

        if self.FREQAXIS == 0:
            for idx in freq_bins:
                self.dynspec[idx,:] = mask_value
                self.mask[idx,:] = True
        else: # FREQAXIS == 1
            for idx in freq_bins:
                self.dynspec[:,idx] = mask_value
                self.mask[:,idx] = True

    def prune_time(self, nsecs, before=True, after=True):
        nbins = int(nsecs/self.dt)
        if before:
            self.dynspec = np.delete(self.dynspec, range(nbins), axis=self.TIMEAXIS)
            self.mask = np.delete(self.mask, range(nbins), axis=self.TIMEAXIS)
            self.t = np.delete(self.t, range(nbins))
            if self.fscrunched is not None:
                self.fscrunched = np.delete(self.fscrunched, range(nbins))
            self.Nt -= nbins
        if after:
            self.dynspec = np.delete(self.dynspec, range(self.Nt-nbins, self.Nt), axis=self.TIMEAXIS)
            self.mask = np.delete(self.mask, range(self.Nt-nbins, self.Nt), axis=self.TIMEAXIS)
            self.t = np.delete(self.t, range(self.Nt-nbins, self.Nt))
            if self.fscrunched is not None:
                self.fscrunched = np.delete(self.fscrunched, range(self.Nt-nbins, self.Nt))
            self.Nt -= nbins

    def add_dm_padding(self, DM, fill_value=None, mask=True):
        before_padding = calc_dmdelay(DM, self.f[0], self.freq_ref)
        after_padding = calc_dmdelay(DM, self.freq_ref, self.f[-1])

        self.add_time_padding(before_padding, fill_value=fill_value, before=True, after=False)
        self.add_time_padding(after_padding, fill_value=fill_value, before=False, after=True)

    def add_time_padding(self, nsecs, fill_value=None, before=True, after=True, mask=True):
        nbins = int(nsecs/self.dt)
        if fill_value is None:
            fill_value = 0.0
        padding_shape = (nbins, self.Nf) if self.TIMEAXIS == 0 else (self.Nf, nbins)
        if before:
            self.dynspec = np.concatenate((np.full(padding_shape, fill_value), self.dynspec), axis=self.TIMEAXIS)
            self.mask = np.concatenate((np.full(padding_shape, mask), self.mask), axis=self.TIMEAXIS)
            if self.fscrunched is not None:
                self.fscrunched = np.concatenate((np.full(nbins, fill_value), self.fscrunched))
            self.Nt += nbins
            self.create_time_axis(self.dt, t0=self.t[0] - self.dt/2 - nbins*self.dt)
        if after:
            self.dynspec = np.concatenate((self.dynspec, np.full(padding_shape, fill_value)), axis=self.TIMEAXIS)
            self.mask = np.concatenate((self.mask, np.full(padding_shape, mask)), axis=self.TIMEAXIS)
            if self.fscrunched is not None:
                self.fscrunched = np.concatenate((self.fscrunched, np.full(nbins, fill_value)))
            self.Nt += nbins
            self.create_time_axis(self.dt, t0=self.t[0] - self.dt/2)

    def transpose(self):
        self.dynspec = self.dynspec.T
        self.mask = self.mask.T

    def create_time_axis(self, sample_time, t0):
        self.t = np.arange(self.Nt)*sample_time + t0 + sample_time/2
        self.dt = sample_time

    def create_freq_axis(self, bw, freqlo):
        self.f = np.arange(self.Nf)*bw + freqlo
        self.df = bw

    def plot(self, ax):
        extent = (
                self.t[0]  - self.dt/2,
                self.t[-1] + self.dt/2,
                self.f[0]  - self.df/2,
                self.f[-1] + self.df/2)
        # Create masked array
        mask_applied = copy.deepcopy(self.dynspec)
        mask_applied[self.mask] = np.nan
        ax.imshow(mask_applied, aspect='auto', interpolation='none', origin='lower', extent=extent)
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

    def calc_dm_shifts(self, dm):
        if self.freq_ref is None:
            raise Exception("Cannot calculate DM shifts without first setting a reference frequency")

        # Calculate how much the DM is to change by
        dm_diff = dm - self.dm

        # The DM delay formula assumes time in seconds, frequencies in MHz
        return calc_dmdelay(dm_diff, self.freq_ref, self.f)

    def dedisperse(self, dm):
        shifts = self.calc_dm_shifts(dm)

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

        # The mask, however, cannot be FFT'd, so we just have to roll each row by an appropriate amount
        for i in range(self.Nf):
            if self.FREQAXIS == 0:
                self.mask[i,:] = np.roll(self.mask[i,:], int(np.round(shifts[i]/self.dt)))
            else: # self.FREQAXIS == 1:
                self.mask[:,i] = np.roll(self.mask[:,i], int(np.round(shifts[i]/self.dt)))

        # Record the current DM
        self.dm = dm

    def fscrunch(self):
        mask_applied = copy.deepcopy(self.dynspec)
        mask_applied[self.mask] = np.nan
        self.fscrunched = np.nanmean(mask_applied, axis=self.FREQAXIS)

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
            self.dynspec.dedisperse(dm)
            self.dynspec.fscrunch()
            self.peak_snrs.append(np.nanmax(self.dynspec.fscrunched))
        self.peak_snrs = np.array(self.peak_snrs)

    def calc_best_dm(self):
        '''
        Assumes run_dmtrials has already been called
        '''
        # Interpolate the DM curve and get the max value
        interped = interp1d(self.dms, self.peak_snrs, kind='cubic')
        #closest_best_dm_idx = np.argmax(self.peak_snrs, keepdims=True)[0]
        self.best_dm = fminbound(lambda dm: -interped(dm), self.dms[0], self.dms[-1])

def main(**kwargs):
    # Load the data into a "dynspec" object
    dynspec = Dynspec(**kwargs)

    # Set the reference frequency
    dynspec.set_freq_ref(kwargs['freq_ref'])

    # Apply padding
    if kwargs['padding'] is not None:
        if kwargs['padding'] == "DM":
            if len(kwargs['dms']) == 1:
                DM = kwargs['dms'][0]
            else:
                DM = kwargs['dms'][1]
            dynspec.add_dm_padding(DM, fill_value=kwargs['mask_value'])
        else:
            try:
                dynspec.add_time_padding(float(kwargs['padding']), fill_value=kwargs['mask_value'])
            except:
                raise Exception("Could not interpret {} as a valid padding amount".format(kwargs['padding']))

    if len(kwargs['dms']) == 1:
        # In this case, only a single DM was given, so we don't
        # do any search, and just use this DM for dedispersion
        DM = kwargs['dms'][0]
    else:
        # Run a bunch of DM trials to create a DM curve
        dmcurve = DMCurve(dynspec)
        dmcurve.run_dmtrials(kwargs['dms'], freq_ref=kwargs['freq_ref'])
        dmcurve.calc_best_dm()

        if kwargs['dmcurve_image'] is not None:
            # Plot the DM curve
            fig, ax = plt.subplots(nrows=1, ncols=1)
            ax.plot(dmcurve.dms, dmcurve.peak_snrs)
            ax.set_xlabel("DM (pc/cm^3)")
            ax.set_ylabel("Peak flux density (a.u.)")

            if kwargs['dmcurve_image'] == "SHOW":
                plt.show()
            else:
                plt.savefig(kwargs['dmcurve_image'])

        DM = dmcurve.best_dm

    # Dedisperse the spectrum
    dynspec.dedisperse(DM)
    dynspec.fscrunch()

    # Plot the dynamic spectrum at the given/best DM
    if kwargs['dynspec_image'] is not None:
        dynspec.t = dynspec.get_time_at_infinite_frequency()
        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
        dynspec.plot_lightcurve(axs[0])
        dynspec.plot(axs[1])
        fig.suptitle('DM = {:.1f} pc/cm^3'.format(DM))
        axs[0].set_yticks([])
        if kwargs['dynspec_image'] == "SHOW":
            plt.show()
        else:
            plt.savefig(kwargs['dynspec_image'])

    # If requested, write out a time series of the frequency-scrunched
    # lightcurve
    if kwargs['lightcurve'] is not None:
        header = "Created with:\n"
        header += '  ' + ' '.join(sys.argv) + '\n\n'
        header += 'This time series has been dedispersed to {} pc/cm^3\n'.format(DM)
        # Get the time of the first bin referenced to infinite frequency
        timeaxis = dynspec.get_time_at_infinite_frequency()
        if kwargs['bc_corr'] == True:
            if kwargs['RA'] is not None and kwargs['Dec'] is not None:
                from bc_corr import bc_corr
                coord = SkyCoord(ra=kwargs['RA']*u.hr, dec=kwargs['Dec']*u.deg, frame='icrs')
                time = Time(timeaxis[0], format='gps')
                bc_correction = bc_corr(coord, time, EPHEMERIS)
                header += 'Using barycentric correction of {} s\n'.format(bc_correction)
                timeaxis += bc_correction # Add the barycentric correction
            else:
                raise Exception("Barycentric correction requested but source coordinates not given")

        # Create verbose header for lightcurve output files
        header += 'Using dedispersion delay of {} s (for reference frequency {} )\n\n'.format(dynspec.dmdelay, dynspec.freq_ref)
        header += "Time (s) | Flux density (a.u.)"

        # Construct array to be written out and write it out
        lightcurve = np.array([timeaxis, dynspec.fscrunched]).T
        np.savetxt(kwargs['lightcurve'], lightcurve, header=header)

    if kwargs['output'] is not None:
        np.savetxt(kwargs['output'], dynspec.dynspec)

def parse_yaml(yaml_file):
    '''
    yaml_file = open file stream for yaml file
    Returns dictionary with abbreviated keys
    '''
    obj = {}

    yaml_params = yaml.safe_load(yaml_file)
    try:
        obj['bc_corr'] = yaml_params['Apply barycentric correction']
    except:
        pass
    try:
        obj['freqlo'] = yaml_params['Dynamic spectrum']['Centre of lowest channel (MHz)']
    except:
        pass
    try:
        obj['bw'] = yaml_params['Dynamic spectrum']['Channel width (MHz)']
    except:
        pass
    try:
        obj['input'] = yaml_params['Dynamic spectrum']['Input file']
    except:
        pass
    try:
        obj['sample_time'] = yaml_params['Dynamic spectrum']['Sample time (s)']
    except:
        pass
    try:
        obj['t0'] = yaml_params['Dynamic spectrum']['T0 (s)']
    except:
        pass
    try:
        obj['transpose'] = yaml_params['Dynamic spectrum']['Transpose']
    except:
        pass
    try:
        obj['freq_ref'] = yaml_params['Reference frequency (MHz)']
    except:
        pass
    try:
        obj['mask_time_bins'] = yaml_params['RFI Mask']['Time bins']
    except:
        pass
    try:
        obj['mask_freq_bins'] = yaml_params['RFI Mask']['Freq bins']
    except:
        pass
    try:
        obj['mask_value'] = yaml_params['RFI Mask']['Value']
    except:
        pass
    try:
        obj['RA'] = yaml_params['RA']
    except:
        pass
    try:
        obj['Dec'] = yaml_params['Dec']
    except:
        pass
    try:
        obj['padding'] = yaml_params['Padding']
    except:
        pass

    return obj

if __name__ == "__main__":
    # Parse the command line
    parser = argparse.ArgumentParser(description='Dedisperse a dynamic spectrum')
    parser.add_argument('--dms', type=float, nargs='*', help='DM trials (pc/cm^3) [[start=0], stop, [step=1]] (i.e. same argument structure as s NumPy''s arange function). If a single value is given, the dynamic spectrum is dedispersed to this DM and no search is performed')
    parser.add_argument('--sample_time', type=float, help='The time of one sample (s)')
    parser.add_argument('--freqlo', type=float, help='The centre frequency of the lowest channel (MHz)')
    parser.add_argument('--freq_ref', type=str, help='The reference frequency used during dedispersion. Either a frequency in MHz, or one of [\'low\', \'centre\', \'high\']')
    parser.add_argument('--bw', type=float, help='The channel width (MHz)')
    parser.add_argument('--output', type=argparse.FileType('w'), help='The file to which the dedispersed dynamic spectrum data will be written. If set to SHOW, then the image is shown (plt.show()) instead.')
    parser.add_argument('--dynspec_image', type=str, help='The file to which the dedispersed dynamic spectrum image will be saved. If set to SHOW, then the image is shown (plt.show()) instead.')
    parser.add_argument('--dmcurve_image', type=str, help='The file to which the DM curve image will be saved')
    parser.add_argument('--input', type=str, help='The (NumPy-readable) file containing the input dynamic spectrum')
    parser.add_argument('--transpose', help='Interpret the input file as rows for time axis, columns for frequency axis')
    parser.add_argument('--lightcurve', type=argparse.FileType('w'), help='Write out the frequency-scrunched, dispersion-corrected lightcurve to the named file. "Dispersion-corrected" means using infinite frequency as reference')
    parser.add_argument('--t0', type=float, help='The left (early) edge of the first time bin')
    parser.add_argument('--bc_corr', help='Apply barycentric correction')
    parser.add_argument('--mask_time_bins', type=int, nargs='*', help='Mask these time bins (expecting ints)')
    parser.add_argument('--mask_freq_bins', type=int, nargs='*', help='Mask these frequency bins (expecting ints)')
    parser.add_argument('--mask_value', type=float, help='The value to use for masked bins/frequencies')
    parser.add_argument('--padding', type=str, help='Number of seconds worth of padding to use along the time axis before dedispersion. The padded pixels will be filled with the value of --mask_value. If set to "DM", and if only one DM is given, then the padding will be chosen to just ensure that there is no wrapping during dedispersion')
    parser.add_argument('--RA', type=float, help='The RA of the source in decimal hours')
    parser.add_argument('--Dec', type=float, help='The Dec of the source in decimal degrees')
    parser.add_argument('--yaml', type=argparse.FileType('r'), help='Obtain parameters from yaml file. These will be overriden by equivalent parameters given on the command line')
    parser.add_argument('--yaml_help', action='store_true', help='More detailed documentation on the --yaml option')

    args = parser.parse_args()

    if args.yaml_help == True:
        print("Documentation for the --yaml option")
        print("===================================\n")
        print("yaml files can be used in place of command line options as an effective 'configuration' file for processing")
        print("dynamic spectra. Each command line option maps to a yaml field, as detailed below. Note that some fields must")
        print("be in particular heirarchies, i.e. as subfields of other fields.\n")
        print("yaml field                     | parent field     | valid values   | equivalent command line option")
        print("-------------------------------+------------------+----------------+-------------------------------")
        print("Apply barycentric correction   |                  | true/false     | --bc_corr")
        print("Dynamic spectrum               |                  |                |")
        print("Centre of lowest channel (MHz) | Dynamic spectrum | [float]        | --freqlo")
        print("Channel width (MHz)            | Dynamic spectrum | [float]        | --bw")
        print("Input file                     | Dynamic spectrum | [file]         | --input")
        print("Sample time (s)                | Dynamic spectrum | [float]        | --sample_time")
        print("T0 (s)                         | Dynamic spectrum | [float]        | --t0")
        print("Transpose                      | Dynamic spectrum | true/false     | --transpose")
        print("Reference frequency (MHz)      |                  | [float]        | --freq_ref")
        print("RA                             |                  | [float]        | --RA")
        print("Dec                            |                  | [float]        | --Dec")
        print("RFI Mask                       |                  |                |")
        print("Time bins                      | RFI Mask         | [list of ints] | --mask_time_bins")
        print("Freq bins                      | RFI Mask         | [list of ints] | --mask_freq_bins")
        print("Padding                        |                  | [float]/DM     | --padding")
        print("\nExample:\n")
        print("Apply barycentric correction: true")
        print("Dynamic spectrum:")
        print("  Centre of lowest channel (MHz): 103.115")
        print("  Channel width (MHz): 0.16")
        print("  Input file: 1343567456_dyn_dynamic_spectrum.csv")
        print("  Sample time (s): 0.5")
        print("  T0 (s): 1343567456.0")
        print("  Transpose: true")
        print("ObsID: 1343567456")
        print("Reference frequency (MHz): centre")
        print("Telescope: MWA")
        print("RFI Mask:")
        print("  Time bins:")
        print("    - 14")
        print("    - 15")
        print("    - 28")
        print("RA: 12.345")
        print("Dec: -67.890")

        exit()

    if args.yaml is not None:
        params = parse_yaml(args.yaml)
        for k, v in vars(args).items():
            if v is not None or k not in params:
                params[k] = v

    main(**params)

