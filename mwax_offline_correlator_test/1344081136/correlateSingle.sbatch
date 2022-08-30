#!/bin/bash --login

#SBATCH --account=mwavcs
#SBATCH --job-name=offlinecorrelate
#SBATCH --mem=370gb
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --partition=gpuq

NUMTILES=$1
MODULE=$(echo "(${NUMTILES} + 8)/16*16" | bc)T

module use /pawsey/mwa/software/python3/modulefiles
module load mwax_offline_correlator/${MODULE}

CHAN=$2
PADDED_CHAN=$(seq -w ${CHAN} 999 | head -1)
srcPath=$3
destinationPath=$4

INPUT_BUFFER_SIZE=$(echo "$NUMTILES * 256000" | bc)
OUTPUT_BUFFER_SIZE=$(echo "($NUMTILES + 1) * $NUMTILES * 102408" | bc)

INPUT_KEY="1${PADDED_CHAN}${PADDED_CHAN}1"
OUTPUT_KEY="2${PADDED_CHAN}${PADDED_CHAN}2"

#for each line/channels make a new buffer
dada_db -b ${INPUT_BUFFER_SIZE} -k ${INPUT_KEY} -n 644 -l -p
dada_db -b ${OUTPUT_BUFFER_SIZE} -k ${OUTPUT_KEY} -n 64 -l -p

#once all started create set proccess for each
mwax_db2fits -k ${OUTPUT_KEY} -l 0 --destination-path=${destinationPath} --health-netiface=lo --health-ip=224.0.2.2 --health-port=8005 > ${destinationPath}/db2fits${PADDED_CHAN}.log 2>&1 &
mwax_db2correlate2db ${INPUT_KEY} ${OUTPUT_KEY} 224.0.2.2 8004 -a 5 -b 160 -d 0 -f 6400 -o 1280 -O 2 -r -v -v > ${destinationPath}/db2c2db${PADDED_CHAN}.log 2>&1 &

cd ${srcPath}
ls *_${CHAN}.sub | sort | sed 's/\(.*\)/dada_diskdb -k '${INPUT_KEY}' -f \1/g' > ch${PADDED_CHAN}.sh
source ch${PADDED_CHAN}.sh > ${destinationPath}/dada_diskdb${PADDED_CHAN}.log 2>&1 &


while [ ! -f ${destinationPath}/*ch${PADDED_CHAN}_000.fits ]; do sleep 1; done
echo "Output FITS file for channel ${CHAN} completed. De-allocating buffer memory..."
dada_db -d -k ${INPUT_KEY}
dada_db -d -k ${OUTPUT_KEY}
