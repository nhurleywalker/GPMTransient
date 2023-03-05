import numpy as np
import matplotlib.pyplot as plt
import dedisperse_dynspec as dd
import pygedm
import astropy.units as u
import copy

from scipy.stats import linregress

# Get the dynamic spectra for the two MWA bands (the "split pulse")
yamlhi = "1342623496.yaml"
with open(yamlhi, 'r') as yaml:
    params = dd.parse_yaml(yaml)
    dhi = dd.Dynspec(**params)

yamllo = "1342623792.yaml"
with open(yamllo, 'r') as yaml:
    params = dd.parse_yaml(yaml)
    dlo = dd.Dynspec(**params)

dhi.set_freq_ref('high')

# Chop off the front chunk of high-frequency data, which has no signal in it
dhi.prune_time(100, before=True, after=False)

# Chop off the tail end of low-frequency data, which has no (or negligible) signal in it
dlo.prune_time(60, before=False, after=True)

# Add enough padding for largest requested DM
dhi.add_dm_padding(350, fill_value=0.0)

# For each candidate DM, dedisperse the high frequency pulse to get the "true pulse shape"
dedispersed_pulse = {}
DMs = np.arange(220, 320)
correlated_NE2001 = []
correlated_YMW16 = []

# In the DM range 250 - 320,
# NE2001 in this direction gives tau_sc = 3.922s x (DM/300)^(3.972) x (f/150 MHz)^(-4.4)
# YMW16  "   "       "       "   tau_sc = 4.005s x (DM/300)^(4.187) x (f/150 MHz)^(-4.0)

for DM in DMs:

    print("Analysing DM = {}...".format(DM))
    # Make a copy and dedisperse it
    dhi_copy = copy.deepcopy(dhi)
    dhi_copy.dedisperse(DM)
    dhi_copy.t = dhi_copy.get_time_at_infinite_frequency()

    # Prune the dedispersed "tails" from the copy
    dhi_copy.prune_time(10, before=True, after=True) # 10 seconds is enough to cover the tails at the highest DM

    # ...and add back a tail of zeros so that the pulse doesn't wrap around
    # back to the beginning when it gets scattered
    dhi_copy.add_time_padding(150, before=False, after=True, mask=False)

    # Get the lightcurve at the higher frequency
    dhi_copy.fscrunch()

    # Get the scattering timescales for each frequency channel in the low frequency observation
    # according to both electron density models
    tau_sc_NE2001 = 3.922 * (DM/300)**(3.972) * (dlo.f/150)**(-4.4)
    tau_sc_YMW16  = 4.005 * (DM/300)**(4.187) * (dlo.f/150)**(-4.0)

    # Construct normalised scattering kernels
    t = np.arange(dhi_copy.Nt) * dhi_copy.dt
    kernel_NE2001 = np.array([np.exp(-t/tau) / np.sum(np.exp(-t/tau)) for tau in tau_sc_NE2001])
    kernel_YMW16 = np.array([np.exp(-t/tau) / np.sum(np.exp(-t/tau)) for tau in tau_sc_YMW16])

    # In the above kernel matrices, the frequency axis is 0, time axis is 1
    FREQAXIS = 0
    TIMEAXIS = 1

    # Fourier transform the scattering kernels to get the equivalent transfer function
    H_NE2001 = np.fft.rfft(kernel_NE2001, axis=TIMEAXIS)
    H_YMW16 = np.fft.rfft(kernel_YMW16, axis=TIMEAXIS)

    # Fourier transform the high frequency lightcurve to get its spectrum
    fscrunched_spectrum = np.fft.rfft(dhi_copy.fscrunched)

    # Multiply them together, with broadcasting, so that the spectrum is hit with each kernel
    scattered_spectrum_NE2001 = fscrunched_spectrum[np.newaxis,:] * H_NE2001
    scattered_spectrum_YMW16 = fscrunched_spectrum[np.newaxis,:] * H_YMW16

    # And back to the time domain...
    scattered_NE2001 = np.fft.irfft(scattered_spectrum_NE2001, n=dhi_copy.Nt, axis=TIMEAXIS)
    scattered_YMW16 = np.fft.irfft(scattered_spectrum_YMW16, n=dhi_copy.Nt, axis=TIMEAXIS)

    # Now that the pulse has been scattered, we need to disperse it.
    # By far the easiest way is to put the scattered pulses into a Dynspec struct.
    # It's not the most efficient, but there's no point reinventing the wheel
    d_NE2001 = dd.Dynspec(dynspec=scattered_NE2001,
            sample_time=dhi_copy.dt,
            t0=dhi_copy.t[0] - dhi_copy.dt/2,
            bw=dlo.df,
            freqlo=dlo.f[0])

    d_YMW16 = dd.Dynspec(dynspec=scattered_YMW16,
            sample_time=dhi_copy.dt,
            t0=dhi_copy.t[0] - dhi_copy.dt/2,
            bw=dlo.df,
            freqlo=dlo.f[0])

    # Disperse them (i.e. dedisperse them with a negative DM
    d_NE2001.set_freq_ref("high")
    d_YMW16.set_freq_ref("high")
    d_NE2001.add_time_padding(130, before=False, after=True, mask=False) # 130s is enough for the highest DM we're considering here
    d_YMW16.add_time_padding(130, before=False, after=True, mask=False)
    d_NE2001.dedisperse(-DM)
    d_YMW16.dedisperse(-DM)

    # Add the bulk dispersion time shift for the reference frequency
    d_NE2001.t += dd.calc_dmdelay(DM, d_NE2001.freq_ref, np.inf)
    d_YMW16.t += dd.calc_dmdelay(DM, d_YMW16.freq_ref, np.inf)

    # Now, the time bins of the predicted dynamic spectra won't exactly line up
    # with the time bins of the low frequency dynamic spectrum, but I think the
    # (random) error we get by effectively "shifting" the spectra by less than
    # one time bin will be much less significant than the correlated noise.

    # Prune the predicted spectra so that they match the dimensions of the low-frequency
    # spectrum. The following assumes that the sampling time is the same for both
    prune_nsecs_before_NE2001 = np.round((dlo.t[0] - d_NE2001.t[0])/d_NE2001.dt) * d_NE2001.dt
    prune_nsecs_after_NE2001 = np.round((d_NE2001.t[-1] - dlo.t[-1])/d_NE2001.dt) * d_NE2001.dt
    d_NE2001.prune_time(prune_nsecs_before_NE2001, before=True, after=False)
    d_NE2001.prune_time(prune_nsecs_after_NE2001, before=False, after=True)

    prune_nsecs_before_YMW16 = np.round((dlo.t[0] - d_YMW16.t[0])/d_YMW16.dt) * d_YMW16.dt
    prune_nsecs_after_YMW16 = np.round((d_YMW16.t[-1] - dlo.t[-1])/d_YMW16.dt) * d_YMW16.dt
    d_YMW16.prune_time(prune_nsecs_before_YMW16, before=True, after=False)
    d_YMW16.prune_time(prune_nsecs_after_YMW16, before=False, after=True)

    # Now, correlate them!
    correlated_NE2001.append(np.nanmean(d_NE2001.dynspec * dlo.dynspec))
    correlated_YMW16.append(np.nanmean(d_YMW16.dynspec * dlo.dynspec))

#fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True)
#fig.suptitle("DM = {:.2f}".format(DM))

#dhi_copy.plot_lightcurve(ax[0])
#dhi_copy.plot(ax[1])

#plt.xlim([dhi.t[0], dlo.t[-1]])
#plt.savefig('lineup.png')
#dlo.plot(ax[2])
#d_NE2001.plot(ax[0])
#d_YMW16.plot(ax[1])

plt.plot(DMs, correlated_NE2001, label="NE2001")
plt.plot(DMs, correlated_YMW16, label="YMW16")
plt.legend()
plt.xlabel("DM (pc/cm^3)")
plt.ylabel("Correlation (a.u.)")

plt.show()
