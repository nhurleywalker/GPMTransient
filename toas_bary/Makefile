LIGHTCURVES_MOD = $(shell ls ../lightcurves_bary/*_lightcurve.txt)
YAML_FILES = $(patsubst ../lightcurves_bary/%_lightcurve.txt,../lightcurves_nobary/%.yaml,$(LIGHTCURVES_MOD))
TOAS_MOD = $(patsubst ../lightcurves_bary/%_lightcurve.txt,%.tim,$(LIGHTCURVES_MOD))

TOA_ERR_SEC = 300

BEST_DM = $(shell cat ../BEST_DM)

all_toas.tim: $(TOAS_MOD)
	echo "FORMAT 1" > $@
	cat $^ | awk '{$$4 = 38000000; print}' >> $@

#------------------------------------------

%.tim: ../lightcurves_bary/%_lightcurve.txt generate_TOA.py
	python $(word 2,$^) ../lightcurves_nobary/$*.yaml $< $(TOA_ERR_SEC) --toa_precision=0.5 --save_plot=$*_toa.png > $@
