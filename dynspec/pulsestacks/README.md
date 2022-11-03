# Pulsestacks

The PNG is displayed below.

There is also an SVG in this folder.
The annotations on the SVG plot contain relative links to the corresponding dynamic spectra, but they aren't clickable when viewing on GitHub.
But, if you open the downloaded SVG as a local file in the browser, the links should work.

**NOTE**: The links are merely relative paths to the dedispersed spectra images, which must obviously exist first.
They can be created by navigating to the parent directory and running
```
mkdir -p dedispersed_spectra
make all_dedispersed_spectra
```

## Making the pulsestacks (new method with `fold.py`)

Unlike the old method (see below), this method uses $\dot{P}$ as well as $P$ for folding.

To use `fold.py`, you need dedispersed lightcurves and a pulsar-style ephemeris (e.g. [gpm.par](gpm.par)).
The, assuming the (barycentred, dedispersed to infinite frequency) lightcurves are in the parent folder, you would plot the pulsestack with
```
python fold.py --add_phase 0.18 gpm.par ../*_lightcurve.txt
```

## Making the pulsestacks (old method)

In the parent folder, open the `Makefile` and set the `DM` variable to the desired DM.
Then, run
```
touch *.yaml
make -j4 all_lightcurves
```
to remake all the lightcurves, and
```
make pulsestacks/pulsestack_<P>s_DM_<DM>.png
```
where `<P>` is the desired folding period in seconds, and `<DM>` must match the `DM` variable set in `Makefile`.

### Making interactive SVG pulsestacks

(Also in the parent directory to this one:)
```
make pulsestacks/pulsestack_<P>s_DM_<DM>.svg
```

## Dedispersed to 275 pc/cm^3

![1318.2 seconds, DM = 275](pulsestack_1318.2s_DM_275.png)

