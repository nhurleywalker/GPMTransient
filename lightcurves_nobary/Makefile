SRCDIR = ../dynspec
YAML_FILES = $(shell ls $(SRCDIR)/*.yaml | sed 's_^$(SRCDIR)/__')
LIGHTCURVES = $(patsubst %.yaml,%_lightcurve.txt,$(YAML_FILES))
BEST_DM = $(shell cat ../BEST_DM)

EPHEMERIS = ../dynspec/de430.bsp

lightcurves: $(LIGHTCURVES)

yaml_files: $(YAML_FILES)

#------------------------------------------

%.yaml:
	awk -F '[ ]' '/barycentric/ { $$4 = "false" } /Input file/ {$$5 = "$(SRCDIR)/" $$5} {print}' < $(SRCDIR)/$@ > $@

%_lightcurve.txt: ../dynspec/dedisperse_dynspec.py %.yaml $(EPHEMERIS)
	python $< --dms $(BEST_DM) --lightcurve $@ --yaml $(word 2,$^)

