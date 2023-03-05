YAML_FILES = $(shell ls ../*.yaml | sed 's_../__')
GPSTIMES = $(patsubst %.yaml,%,$(YAML_FILES))
LIGHTCURVES = $(patsubst %,%_lightcurve.txt,$(GPSTIMES))
TOAS = $(patsubst %,%.tim,$(GPSTIMES))

TOA_ERR_SEC = 300

BEST_DM = 275

prepare_yaml_files: $(YAML_FILES)

prepare_lightcurves: $(LIGHTCURVES)

all_toas.tim: $(TOAS)
	echo "FORMAT 1" > $@
	cat $^ | awk '{$$4 = 60000000; print}' >> $@

#------------------------------------------

%.yaml: ../%.yaml
	awk -F '[ ]' '/barycentric/ { $$4 = "false" } /Input file/ {$$5 = "../" $$5} {print}' < $< > $@

%_lightcurve.txt: ../dedisperse_dynspec.py %.yaml ../de430.bsp
	python $< --dms $(BEST_DM) --lightcurve $@ --yaml $(word 2,$^)

%.tim: %_lightcurve.txt generate_TOA.py
	python $(word 2,$^) $*.yaml $< $(TOA_ERR_SEC) --toa_precision=0.5 --save_plot=$*_toa.png > $@
