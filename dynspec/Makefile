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
	python $< $$(awk '$$1 == $* {if ($$6 == "Y") {printf "--transpose "}; printf "--sample_time %s --bw %s --freqlo %s --dms 150 401 5 --input %s --lightcurve $@ --t0 $* --no_plots\n", $$3, $$4, $$5, $$7}' < $(word 2,$^))
