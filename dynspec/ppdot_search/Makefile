YAML_FILES = $(shell ls ../*.yaml | sed 's_../__')
GPSTIMES = $(patsubst %.yaml,%,$(YAML_FILES))
LIGHTCURVES = $(patsubst %,%_lightcurve.txt,$(GPSTIMES))

BEST_DM = 275

prepare_yaml_files: $(YAML_FILES)

prepare_lightcurves: $(LIGHTCURVES)

#------------------------------------------

%.yaml: ../%.yaml
	awk -F '[ ]' '/barycentric/ { $$4 = "false" } /Input file/ {$$5 = "../" $$5} {print}' < $< > $@

%_lightcurve.txt: ../dedisperse_dynspec.py %.yaml ../de430.bsp
	python $< --dms $(BEST_DM) --lightcurve $@ --yaml $(word 2,$^)
