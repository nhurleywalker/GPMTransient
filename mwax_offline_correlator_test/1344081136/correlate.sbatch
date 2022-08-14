#!/bin/bash -l

#SBATCH --account=mwavcs
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --output=mwax_offline_correlate_1344081136-%j.out
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:1
#SBATCH --mem=370gb

# At this time of writing, this should NOT be run as a batch job,
# but rather run on an interactive node (via salloc) so that the
# jobs finish cleanly before the memory is cleared up at the end.

NUMTILES=136

module use /pawsey/mwa/software/python3/modulefiles
#module load mwax_offline_correlator/${NUMTILES}T
module load mwax_offline_correlator/144T

INPUT_BUFFER_SIZE=$(echo "$NUMTILES * 256000" | bc)
OUTPUT_BUFFER_SIZE=$(echo "($NUMTILES + 1) * $NUMTILES * 102408" | bc)
dada_db -b $INPUT_BUFFER_SIZE -k 1234 -n 644 -l -p
dada_db -b $OUTPUT_BUFFER_SIZE -k 2345 -n 64 -l -p

mwax_db2fits -k 2345 --destination-path=/astro/mwavcs/smcsweeney/1344081136_correlated --health-netiface=lo --health-ip=224.0.2.2 --health-port=8005 > db2fits.log 2>&1 &
mwax_db2correlate2db 1234 2345 224.0.2.2 8004 -a 5 -b 160 -d 0 -f 6400 -o 1280 -O 2 -r -t -v -v > db2c2db.log 2>&1 &

for CHAN in {109..132}; do
    find /astro/mwavcs/smcsweeney/1344081136/ -type f -name "*_${CHAN}.sub" | sort |
        while read subfile; do
            dada_diskdb -k 1234 -f $subfile
        done
done

dada_db -d -k 1234
dada_db -d -k 2345
