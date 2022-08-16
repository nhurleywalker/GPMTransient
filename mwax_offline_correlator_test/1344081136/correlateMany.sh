#!/bin/bash



# Command line arguements are start channel, end channel, input and output locations

NUMTILES=136

startChan=$1
endChan=$2
indir=$3
outdir=$4


INPUT_BUFFER_SIZE=$(echo "$NUMTILES * 256000" | bc)
OUTPUT_BUFFER_SIZE=$(echo "($NUMTILES + 1) * $NUMTILES * 102408" | bc)



MODULE=$(echo "(${NUMTILES} + 8)/16*16" | bc)T

module use /pawsey/mwa/software/python3/modulefiles
module load mwax_offline_correlator/${MODULE}



#for each line/channels make a new buffer
for CHAN in $(seq $startChan 1 $endChan); do

dada_db -b 34816000 -k "1${CHAN}${CHAN}1" -n 644 -l -p
dada_db -b 1908065856 -k "2${CHAN}${CHAN}2" -n 64 -l -p

done


#once all started create set proccess for each
for CHAN in $(seq $startChan 1 $endChan); do

mwax_db2fits -k "2${CHAN}${CHAN}2" --destination-path=outdir -l 0 --health-netiface=lo --health-ip=224.0.2.2 --health-port=8005 > db2fits${CHAN}.log 2>&1 &
mwax_db2correlate2db "1${CHAN}${CHAN}1" "2${CHAN}${CHAN}2" 224.0.2.2 8004 -a 5 -b 160 -d 0 -f 6400 -o 1280 -O 2 -r -v -v > db2c2db${CHAN}.log 2>&1 &

done



s

cd ${indir}

ls *_${CHAN}.sub | sort | sed 's/\(.*\)/srun dada_diskdb -k 1'${CHAN}${CHAN}'1 -f \1/g' > ch${CHAN}.sh
source ch${CHAN}.sh > ${destinationPath}/dada_diskdb${CHAN}.log 2>&1 &


#kill buffers once the fits file is present
for CHAN in $(seq $startChan 1 $endChan); do
    while [ ! -f ${4}/*ch${CHAN}_000.fits ]; do sleep 1; done
    dada_db -d -k "1${CHAN}${CHAN}1"
    dada_db -d -k "2${CHAN}${CHAN}2"
done 


