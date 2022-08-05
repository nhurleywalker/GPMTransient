GPSTIMES = $(shell grep "^[0-9]" observations.txt | cut -d' ' -f1)
BEST_DM = 285

# Lightcurves DM-corrected to infinite frequency, and also barycentric corrected
LIGHTCURVES = $(patsubst %,%_lightcurve.txt,$(GPSTIMES))

.SECONDARY:

all_lightcurves: $(LIGHTCURVES)

%_lightcurve.txt: dedisperse_dynspec.py observations.txt
	python $< $$(awk '$$1 == $* {t0 = $* + $$9; if ($$6 == "Y") {printf "--transpose "}; printf "--sample_time %s --bw %s --freqlo %s --dms $(BEST_DM) --input %s --bc_corr %s --lightcurve $@ --t0 %f --no_plots\n", $$3, $$4, $$5, $$7, $$8, t0}' < $(word 2,$^))

# This recipe is here just as an example of a call to bc_corr.py.
# It won't work "out of the box" here, because the paths to the
# bc_corr binary and also to the de430.bsp file (ephemeris) will be
# different for each user. Run "python bc_corr.py -h" for usage.

bc_corrections: bc_corr.py
	python $< 18.6505 -10.5304 $$(locate de430.bsp | head -1) $(GPSTIMES)

pulsestacks/pulsestack_%s_DM_$(BEST_DM).png: $(LIGHTCURVES)
	python make_pulsestack.py --png $@ $* $^

