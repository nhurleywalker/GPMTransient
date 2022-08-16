import numpy as np
import matplotlib.pyplot as plt
import argparse

if __name__ == '__main__':
    # Parse the command line
    parser = argparse.ArgumentParser(description='Make a pulsestack from a collection of light curves')
    parser.add_argument('period', type=float, help='The folding period in seconds')
    parser.add_argument('lightcurves', type=str, nargs='*', help='The files containing lightcurve data. Expected to contain two columns: (1) Time (s), (2) Flux density')
    parser.add_argument('--time_offset', type=float, help='An optional offset')
    parser.add_argument('--png', type=str, help='Output PNG file.')
    parser.add_argument('--svg', type=str, help='Output SVG file.')
    parser.add_argument('--show_plot', action='store_true', help='Show interactive plot.')

    args = parser.parse_args()

    # Shorthand:
    P = args.period

    if args.png is not None or args.svg is not None:
        fig = plt.figure(figsize=(8,16))

    for i in range(len(args.lightcurves)):
        file = args.lightcurves[i]
        obsname = file[:10]
        lightcurve = np.loadtxt(file)
        gpstimes = lightcurve[:,0] # First column
        if i == 0:
            yticks = [0]
            pulse_numbers = ["0"]
            gps_ref = gpstimes[0]
        else:
            pulse_number = int(np.round((gpstimes[0] - gps_ref)/P))
            if pulse_number != pulse_numbers[-1]:
                yticks.append(yticks[-1] + 1)
                pulse_numbers.append(pulse_number)

        phase = (gpstimes - gps_ref + P/2) % P - P/2
        flux_density = lightcurve[:,1] # Second column
        x = phase
        y = 0.7*flux_density/np.max(flux_density) + yticks[-1]
        plt.plot(x, y, lw=0.5)
        url = "../dedispersed_spectra/" + obsname + "_dedispersed.png"
        plt.text(x[0], y[0], obsname, url=url, bbox = dict(color='w', alpha=0.01, url=url))
        plt.yticks(ticks=yticks, labels=pulse_numbers)

    plt.xlabel("Time (s)")
    plt.ylabel("Pulse number")
    plt.title("Period = {:.3f} s".format(P))

    if args.png is not None:
        plt.savefig(args.png)
    if args.svg is not None:
        plt.savefig(args.svg)

    if args.show_plot == True:
        plt.show()
