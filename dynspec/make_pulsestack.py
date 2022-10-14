import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
from astropy.time import Time
import argparse

import matplotlib.font_manager
from matplotlib import rc
# Nature requires sans-serif fonts
plt.rcParams.update({
    "text.usetex": False,
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

# inches to cm
cm = 1/2.54

if __name__ == '__main__':
    # Parse the command line
    parser = argparse.ArgumentParser(description='Make a pulsestack from a collection of light curves')
    parser.add_argument('period', type=float, help='The folding period in seconds')
    parser.add_argument('lightcurves', type=str, nargs='*', help='The files containing lightcurve data. Expected to contain two columns: (1) Time (s), (2) Flux density')
    parser.add_argument('--time_offset', type=float, help='An optional offset')
    parser.add_argument('--png', type=str, help='Output PNG file.')
    parser.add_argument('--pdf', type=str, help='Output PDF file.')
    parser.add_argument('--svg', type=str, help='Output SVG file.')
    parser.add_argument('--show_plot', action='store_true', help='Show interactive plot.')
    parser.add_argument('--gpstimes', default=False, action='store_true', help='Show GPS times on the pulse stack.')
    parser.add_argument('--dates', default=False, action='store_true', help='Show dates on the pulse stack (overrides gpstimes)')
    parser.add_argument('--title', default=False, action='store_true', help='Show period in title')

    args = parser.parse_args()

    # Shorthand:
    P = args.period

    if args.png is not None or args.svg is not None or args.pdf is not None:
# https://www.nature.com/nature/for-authors/final-submission#:~:text=For%20guidance%2C%20Nature's%20standard%20figure,(120%E2%80%93136%20mm).
        fig = plt.figure(figsize=(8.9*cm,16*cm))

    cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["c","m"])

# Some pulses have low signal-to-noise or we only catch a very small amount of the pulse; don't display these
    n = 0
    for i in range(len(args.lightcurves)):
        file = args.lightcurves[i]
        obsname = file[:10]
        ti = Time(obsname, format="gps")
        ti.format="ymdhms"
        t = ti.value
        lightcurve = np.loadtxt(file)
        with open(file, 'r') as f:
            for line in f.readlines():
                if 'frequency' in line:
                    freq = float(line.split(' ')[10])
        gpstimes = lightcurve[:,0] # First column
        flux_density = lightcurve[:,1] # Second column
        if np.isnan(flux_density).any():
            print(f"found some nans for {obsname}, not proceeding further")
        elif obsname == "1342096266":
            print("Skipping Parkes data")
        else:
            if i == 0:
                yticks = [0]
                pulse_numbers = ["0"]
                gps_ref = gpstimes[0]
                xt = -580
                t0 = ti
            else:
                pulse_number = int(np.round((gpstimes[0] - gps_ref)/P))
                if pulse_number == pulse_numbers[-1]:
                    xt = 300 
                else:
                    xt = -580
                    yticks.append(yticks[-1] + 1)
                    pulse_numbers.append(pulse_number)
            phase = (gpstimes - gps_ref + P/2) % P - P/2
            x = phase
            y = 0.7*flux_density/np.nanmax(flux_density) + yticks[-1]
# Log shows the frequencies a little better than linear
            plt.plot(x, y[np.isfinite(y)], lw=0.5, color=cm1((np.log10(freq) - np.log10(88.))/(np.log10(500.)-np.log10(88.))))
            #plt.plot(x, y[np.isfinite(y)], lw=0.5, color=cm1((freq - 88.)/(500.-88.)))
            url = "../dedispersed_spectra/" + obsname + "_dedispersed.png"
            n+=1
            if args.dates is True:
                plt.text(xt, yticks[-1], f"{t[0]:04d}-{t[1]:02d}-{t[2]:02d} {t[3]:02d}:{t[4]:02d}", url=url, bbox = dict(color='w', alpha=0.01, url=url), fontsize=5)
            elif args.gpstimes is True:
                plt.text(xt, yticks[-1], obsname, url=url, bbox = dict(color='w', alpha=0.01, url=url))
    # Turn off yticks for now
    #        plt.yticks(ticks=yticks, labels=pulse_numbers, fontsize=5)
    plt.yticks(ticks=[], labels=[], fontsize=5)
    plt.xlim([-600, 600])
    plt.ylim([yticks[0]-1, yticks[-1]+1])

    plt.axvline(-100, lw=0.5, color="k", alpha=0.5, ls="--")
    plt.axvline(250, lw=0.5, color="k", alpha=0.5, ls="--")
    plt.xlabel("Time (s)")
    #plt.ylabel("Pulse number")
    if args.title is True:
        plt.title("Period = {:.5f} s".format(P))

    if args.png is not None:
        plt.savefig(args.png, bbox_inches="tight")
    if args.svg is not None:
        plt.savefig(args.svg)
    if args.pdf is not None:
        plt.savefig(args.pdf, bbox_inches="tight")

    if args.show_plot == True:
        plt.show()

    timedelta = (ti - t0) / 365.25
    print(f"Displayed {n} light curves over {timedelta.value:2.1f} years.")
