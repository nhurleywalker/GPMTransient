import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.font_manager
from matplotlib import rc
# Nature requires sans-serif fonts
plt.rcParams.update({
    "text.usetex": False,
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

# inches to cm
cm = 1/2.54

from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
import astropy.units as u
import argparse

class Ephemeris:

    def __init__(self, parfileobject):
        for line in parfileobject:
            parsed = line.split()

            # Make sure the parsed line always has at least three values
            if len(parsed) < 2:
                continue

            if len(parsed) == 2:
                parsed.append(None)

            setattr(self, parsed[0], {'value': parsed[1], 'error': parsed[2]})

        # Create fields with units (using lower case variable names)
        self.pos = SkyCoord(f'{self.RAJ["value"]} {self.DECJ["value"]}', unit=(u.hourangle, u.deg), frame='icrs')
        self.dm = float(self.DM['value']) * u.parsec / u.cm**3
        self.pepoch = Time(self.PEPOCH['value'], format='mjd')
        self.f0 = float(self.F0['value']) * u.Hz
        self.f1 = float(self.F1['value']) * u.Hz/u.second
        self.posepoch = Time(self.POSEPOCH['value'], format='mjd')
        self.dmepoch = Time(self.DMEPOCH['value'], format='mjd')
        self.p0 = (1/self.f0).decompose()
        self.p1 = (-self.f1/self.f0**2).decompose()

if __name__ == '__main__':
    # Parse the command line
    parser = argparse.ArgumentParser(description='Convert a lightcurve''s time to rotation phase according to a given ephemeris')
    parser.add_argument('par', type=argparse.FileType('r'), help='The ephemeris file to use for folding')
    parser.add_argument('lightcurves', type=str, nargs='*', help='The ASCII files containing lightcurve data. Expected to contain two columns: (1) Time (GPS seconds), (2) Flux density')
    parser.add_argument('--p0_gps', type=float, help='Number the pulses from this GPS time (default: set earliest pulse to pulse #0)')
    parser.add_argument('--add_phase', type=float, default=0, help='Add phase rotation (0 to 1) to each pulse on stack for display (default: 0)')
    parser.add_argument('--plot_image', type=str, help='Write the plot image to the named file. If not provided, the plot will be shown on screen.')

    args = parser.parse_args()

    eph = Ephemeris(args.par)

    # Open the lightcurves
    phs = []
    pulses = []
    Is = []
    freqs = []

    # Collect all the data and calculate the phases, pulse numbers, etc
    npulses = len(args.lightcurves)
    for i in range(npulses):
        data = np.loadtxt(args.lightcurves[i])
        with open(args.lightcurves[i], 'r') as f:
            for line in f.readlines():
                if 'frequency' in line:
                    freqs.append(float(line.split(' ')[10]))
        t = TimeDelta(data[:,0] - eph.pepoch.gps, format='sec').to(u.second)
        I = data[:,1] * u.jansky
        N = eph.f0*t + 0.5*eph.f1*t**2 + args.add_phase
        pulse, ph = np.divmod(N, 1)
        ph -= 0.5

        Is.append(I)
        phs.append(ph)
        pulses.append(int(pulse[len(t)//2]))

    # Sort the pulses by pulse number
    sidxs = np.argsort(pulses)

    # Change pulse numbers relative to the chosen reference
    fig = plt.figure(figsize=(8.9*cm,16*cm))
    cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["c","m"])
    if args.p0_gps:
        t0 = TimeDelta(args.p0_gps - eph.pepoch.gps, format='sec').to(u.second)
        pulse0 = int(np.round(eph.f0*t0 + 0.5*eph.f1*t0**2))
    else:
        pulse0 = pulses[sidxs[0]]

    pulses = [pulse - pulse0 for pulse in pulses]

    # Plot everything
    for i in range(npulses):
        ph = phs[sidxs[i]]
        pulse = pulses[sidxs[i]]
        I = Is[sidxs[i]]
        freq = freqs[sidxs[i]]
        color = cm1((np.log10(freq) - np.log10(88.))/(np.log10(500.)-np.log10(88.)))
        plt.plot(ph, I/np.max(I) + i, lw=0.5, color=color)
    plt.yticks(ticks=range(npulses), labels=pulses)

    if args.plot_image:
        plt.savefig(args.plot_image)
    else:
        plt.show()
