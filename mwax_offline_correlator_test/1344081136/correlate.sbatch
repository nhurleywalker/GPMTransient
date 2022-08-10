#!/bin/bash -l

#SBATCH --account=mwavcs
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --ntasks=3
#SBATCH --output=correlate.out
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:1
#SBATCH --mem=370gb

module use /pawsey/mwa/software/python3/modulefiles
module load mwax_offline_correlator

# input buffer size = NUMTILES * 256000 = 136 * 256000 = 34816000
dada_db -b 34816000 -k 1234 -n 644 -l -p
# output buffer size = (NUMTILES + 1) * NUMTILES * 102408 = 1908065856
dada_db -b 1908065856 -k 2345 -n 64 -l -p

mwax_db2fits -k 2345 --destination-path=/astro/mwavcs/smcsweeney/mwax_offline_correlator_test --health-netiface=lo --health-ip=224.0.2.2 --health-port=8005 &
mwax_db2correlate2db 1234 2345 224.0.2.2 8004 -a 5 -b 160 -d 0 -f 6400 -o 6400 -O 2 -r -v -v &

srun -n 1 dada_diskdb -k 1234 -f /astro/mwavcs/smcsweeney/mwax_offline_correlator_test/1329828304_1329828464_103.sub -s

#dada_db -d -k 1234
#dada_db -d -k 2345
