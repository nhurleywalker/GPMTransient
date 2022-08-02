GPSTIMES = 1340639536 \
		   1340641016 \
		   1341930288 \
		   1341932656 \
		   1341932952 \
		   1341937984 \
		   1341938280 \
		   1341939168 \
		   1341939760 \
		   1342095399 \
		   1342096104 \
		   1342096400 \
		   1342097584

# Lightcurves DM-corrected to infinite frequency, and also barycentric corrected
LIGHTCURVES = $(patsubst %,%_lightcurve.txt,$(GPSTIMES))

all_lightcurves: $(LIGHTCURVES)

%_lightcurve.txt: dedisperse_dynspec.py observations.txt
	python $< $$(awk '$$1 == $* {if ($$6 == "Y") {printf "--transpose "}; printf "--sample_time %s --bw %s --freqlo %s --dms 150 401 5 --input %s --bc_corr %s --lightcurve $@ --t0 $* --no_plots\n", $$3, $$4, $$5, $$7, $$8}' < $(word 2,$^))

# This recipe is here just as an example of a call to bc_corr.py.
# It won't work "out of the box" here, because the paths to the
# bc_corr binary and also to the de430.bsp file (ephemeris) will be
# different for each user. Run "python bc_corr.py -h" for usage.
#
bc_corrections: bc_corr.py
	python $< 18.6505 -10.5304 $$(locate de430.bsp | head -1) $(GPSTIMES)
