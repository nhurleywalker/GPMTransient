YAML_FILES = $(shell ls ../*.yaml | sed 's_../__')

prepare_yaml_files: $(YAML_FILES)

%.yaml: ../%.yaml
	awk -F '[ ]' '/barycentric/ { $$4 = "false" } /Input file/ {$$5 = "../" $$5} {print}' < $< > $@
