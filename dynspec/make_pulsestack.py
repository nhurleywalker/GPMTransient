import numpy as np
import matplotlib.pyplot as plt
import argparse

if __name__ == '__main__':
    # Parse the command line
    parser = argparse.ArgumentParser(description='Make a pulsestack from a collection of light curves')
    parser.add_argument('period', type=float, help='The folding period in seconds')
    parser.add_argument('lightcurves', type=str, nargs='*', help='The files containing lightcurve data. Expected to contain two columns: (1) Time (s), (2) Flux density')
    parser.add_argument('--time_offset', type=float, help='An optional offset')
    parser.add_argument('--png', type=str, help='Output PNG file. If not set, then the plot is just shown, not saved.')

    args = parser.parse_args()

    # Shorthand:
    P = args.period

    if args.png is not None:
        fig = plt.figure(figsize=(8,12))

    for i in range(len(args.lightcurves)):
        file = args.lightcurves[i]
        lightcurve = np.loadtxt(file)
        gpstimes = lightcurve[:,0] # First column
        if i == 0:
            yticks = [0]
            pulsenumber = ["0"]
            gps_ref = gpstimes[0]
        else:
            yticks.append(i)
            pulsenumber.append(int(np.round((gpstimes[0] - gps_ref)/P)))

        phase = (gpstimes - gps_ref + P/2) % P - P/2
        flux_density = lightcurve[:,1] # Second column
        x = phase
        y = flux_density/2 + i
        plt.plot(x, y)
        plt.yticks(ticks=yticks, labels=pulsenumber)

    plt.xlabel("Time (s)")
    plt.ylabel("Pulse number")

    if args.png is not None:
        plt.savefig(args.png)
    else:
        plt.show()
