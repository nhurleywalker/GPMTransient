GPSTIMES = $(shell ls *.yaml | sed 's/\.yaml//')
BEST_DM = 288

# Lightcurves DM-corrected to infinite frequency, and also barycentric corrected
LIGHTCURVES = $(patsubst %,%_lightcurve.txt,$(GPSTIMES))

# Dedispersed spectra
DEDISPERSED_SPECTRA = $(patsubst %,dedispersed_spectra/%_dedispersed.png,$(GPSTIMES))

.SECONDARY:

all_lightcurves: $(LIGHTCURVES)

%_lightcurve.txt: dedisperse_dynspec.py %.yaml
	python $< --dms $(BEST_DM) --lightcurve $@ --yaml $(word 2,$^)

all_dedispersed_spectra: $(DEDISPERSED_SPECTRA)

dedispersed_spectra/%_dedispersed.png: dedisperse_dynspec.py %.yaml
	python $< --dms $(BEST_DM) --dynspec_image $@ --yaml $(word 2,$^)

# This recipe is here just as an example of a call to bc_corr.py.
# It won't work "out of the box" here, because the paths to the
# bc_corr binary and also to the de430.bsp file (ephemeris) will be
# different for each user. Run "python bc_corr.py -h" for usage.

bc_corrections: bc_corr.py
	python $< --ephemeris $$(locate de430.bsp | head -1) 18.6505 -10.5304 $(GPSTIMES)

pulsestacks/pulsestack_%s_DM_$(BEST_DM).png: $(LIGHTCURVES)
	python make_pulsestack.py --png $@ $* $^

