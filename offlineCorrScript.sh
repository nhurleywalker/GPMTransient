#!/bin/bash

module use /pawsey/mwa/software/python3/modulefiles
module load mwax_offline_correlator/144T


# Command line arguements are input and output locations

#/astro/mwavcs/smcsweeney/1344081136/
#/astro/mwavcs/jmoseley/CorrelatorStuff/CorrelatorOut


#parse file count lines/channels

startChan=109
endChan=132

#for each line/channels make a new buffer
for CHAN in $(seq $startChan 1 $endChan); do

dada_db -b 34816000 -k "${CHAN}A" -n 644 -l -p
dada_db -b 1908065856 -k "${CHAN}B" -n 64 -l -p

done


#once all started create set proccess for each
for CHAN in $(seq $startChan 1 $endChan); do

mwax_db2fits -k "${CHAN}B" --destination-path=${2} --health-netiface=lo --health-ip=224.0.2.2 --health-port=8005 > db2fits${CHAN}.log 2>&1 &
mwax_db2correlate2db "${CHAN}A" "${CHAN}B" 224.0.2.2 8004 -a 5 -b 160 -d 0 -f 6400 -o 6400 -O 2 -r -v -v > db2c2db${CHAN}.log 2>&1 &

done

#sleep here maybe?

#for each Channel/lines

cd ${1}
for CHAN in $(seq $startChan 1 $endChan); do 
 $(ls *_${CHAN}.sub | sed 's/\(.*\)/dada_diskdb -k '${CHAN}'A -f \1 \&\&/g' | tr '\n' ' ' | sed 's/\&\s$//');
done


#cd /astro/mwavcs/jmoseley/mwax_offline_correlator_test

for CHAN in $(seq $startChan 1 $endChan); do
    while [ ! -f ${2}/*ch${CHAN}_000.fits ]; do sleep 1; done
    dada_db -d -k "${CHAN}A"
    dada_db -d -k "${CHAN}B"
done 



#run line code that starts things asyncronisly
#srun -n 1 dada_diskdb -k 1234 -f /astro/mwavcs/smcsweeney/mwax_offline_correlator_test/1329828304_1329828464_103.sub -s



#wait till last fits file exits


#end

