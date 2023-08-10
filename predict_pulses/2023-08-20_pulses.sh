#!/bin/bash

python predict_toas.py dofit.par ../toas_nobary/all_toas_mod.tim 60176 60178  --refit --manual_offset 370 --outparfile out.par

# 2023-08-20_pulses.txt
# The above command outputs "out.tim". The file 2023-08-20_pulses.txt was produced by manually selecting just the subset of MJDs in out.tim which occur between 2023-08-20 14:00:00 and 2023-08-20 20:00:00.

# 2023-08-20_pulses.png
cat ../toas_nobary/all_toas_mod.tim out.tim > all.tim
pintk out.par all.tim
# ...and then grabbing a screenshot of the resulting plot.


# After subtracting 2.5 minutes to signal the beginning of the pulse window, the times can be listed in different formats via:
./timeconv $(grep "^[0-9]" 2023-08-20_pulses.txt | cut -f1 -d' ') --to_fmt gps
