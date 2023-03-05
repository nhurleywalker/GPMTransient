export PYTHONPATH=

# Predict when we would see it in August 2022 with the MWA
python3 predict_pulses.py  --starttime 1343220000 --stoptime 1343500000 --obslength 296 --output GPMJ1839_2022-08_MWA.txt --separation 10

# Predict when we would see it inside the list of observations Tracy sent
python3 predict_pulses.py --input=VLITE_list.txt --telescope=VLA --output=VLITE_predictions.txt
