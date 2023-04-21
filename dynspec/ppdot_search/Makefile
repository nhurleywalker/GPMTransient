LIGHTCURVES = $(shell ls ../*_lightcurve.txt | sed 's_../__')
LIGHTCURVES_MOD = $(shell ls ../*_lightcurve_mod.txt | sed 's_../__')
YAML_FILES = $(patsubst %_lightcurve.txt,%.yaml,$(LIGHTCURVES))
TOAS_ORIG = $(patsubst %_lightcurve.txt,%_orig.tim,$(LIGHTCURVES))
TOAS_MOD = $(patsubst %_lightcurve_mod.txt,%_mod.tim,$(LIGHTCURVES_MOD))

TOA_ERR_SEC = 300

BEST_DM = 275

prepare_yaml_files: $(YAML_FILES)

prepare_lightcurves: $(LIGHTCURVES)

all_toas_orig: $(TOAS_ORIG)

all_toas_mod.tim: $(TOAS_MOD)
	echo "FORMAT 1" > $@
	cat $^ | awk '{$$4 = 60000000; print}' >> $@

#------------------------------------------

%.yaml: ../%.yaml
	awk -F '[ ]' '/barycentric/ { $$4 = "false" } /Input file/ {$$5 = "../" $$5} {print}' < $< > $@

%_lightcurve_mod.txt: ../%_lightcurve_mod.txt
	ln -s $< $@

%_lightcurve.txt: ../dedisperse_dynspec.py %.yaml ../de430.bsp
	python $< --dms $(BEST_DM) --lightcurve $@ --yaml $(word 2,$^)

%_mod.tim: %_lightcurve_mod.txt generate_TOA.py
	python $(word 2,$^) $*.yaml $< $(TOA_ERR_SEC) --toa_precision=0.5 --save_plot=$*_toa.png > $@

%_orig.tim: %_lightcurve.txt generate_TOA.py
	python $(word 2,$^) $*.yaml $< $(TOA_ERR_SEC) --toa_precision=0.5 --save_plot=$*_toa.png > $@
