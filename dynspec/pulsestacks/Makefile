IMAGES = orig.png \
		 fitted.png \
		 second_island.png

all_images: $(IMAGES)

%.png: %.par fold.py
	python $(word 2,$^) --plot_image $@ --add_phase 0.18 $< ../*_lightcurve.txt
