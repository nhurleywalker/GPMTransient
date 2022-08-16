#!/bin/bash --login

#SBATCH --account=mwavcs
#SBATCH --job-name=offlinecorrelate
#SBATCH --mem=370gb
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --partition=gpuq

module use /pawsey/mwa/software/python3/modulefiles
module load mwax_offline_correlator/144T

CHAN=$1
srcPath=$2
destinationPath=$3

#for each line/channels make a new buffer
dada_db -b 34816000 -k "1${CHAN}${CHAN}1" -n 644 -l -p
dada_db -b 1908065856 -k "2${CHAN}${CHAN}2" -n 64 -l -p

#once all started create set proccess for each
mwax_db2fits -k "2${CHAN}${CHAN}2" -l 0 --destination-path=${destinationPath} --health-netiface=lo --health-ip=224.0.2.2 --health-port=8005 > ${destinationPath}/db2fits${CHAN}.log 2>&1 &
mwax_db2correlate2db "1${CHAN}${CHAN}1" "2${CHAN}${CHAN}2" 224.0.2.2 8004 -a 5 -b 160 -d 0 -f 6400 -o 1280 -O 2 -r -v -v > ${destinationPath}/db2c2db${CHAN}.log 2>&1 &

cd ${srcPath}
ls *_${CHAN}.sub | sort | sed 's/\(.*\)/dada_diskdb -k 1'${CHAN}${CHAN}'1 -f \1/g' > ch${CHAN}.sh
source ch${CHAN}.sh > ${destinationPath}/dada_diskdb${CHAN}.log 2>&1 &


while [ ! -f ${destinationPath}/*ch${CHAN}_000.fits ]; do sleep 1; done
dada_db -d -k "1${CHAN}${CHAN}1"
dada_db -d -k "2${CHAN}${CHAN}2"

