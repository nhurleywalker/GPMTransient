LIGHTCURVES_MOD = $(shell ls ../lag_analysis/*_lightcurve_mod.txt)
YAML_FILES = $(patsubst ../lag_analysis/%_lightcurve_mod.txt,../lightcurves_nobary/%.yaml,$(LIGHTCURVES_MOD))
TOAS_MOD = $(patsubst ../lag_analysis/%_lightcurve_mod.txt,%_mod.tim,$(LIGHTCURVES_MOD))

TOA_ERR_SEC = 38
TOA_TRUNC_ERR_SEC = 76

BEST_DM = $(shell cat ../BEST_DM)

all_toas_mod.tim: $(TOAS_MOD)
	echo "FORMAT 1" > $@
	cat $^ >> $@

#------------------------------------------

%_mod.tim: ../lag_analysis/%_lightcurve_mod.txt ../scripts/generate_TOA.py
	python $(word 2,$^) ../lightcurves_nobary/$*.yaml $< $(TOA_ERR_SEC) --truncated_error $(TOA_TRUNC_ERR_SEC) --toa_precision=0.5 --save_plot=$*_toa.png > $@

clean:
	$(RM) *.tim *.png
