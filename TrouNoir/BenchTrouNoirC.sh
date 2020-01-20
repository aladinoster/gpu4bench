#!/bin/bash
#SEQ="6 1 14"
SEQ="6 1 13"
HOST=$(hostname)
DATE=$(date "+%Y%m%d-%H%M")
LOGFILE=trou_noir_FP32-${HOST}-${DATE}.log
echo >$LOGFILE
RAIE=CORPS_NOIR
echo Output stored in $LOGFILE
echo -e "Experience : $RAIE" >>$LOGFILE
echo -e "Experience : $RAIE"
seq $SEQ | while read POWER ; do
    SIZE=$((2**$POWER))
    seq 1 1 10 | xargs -I TOTO  /usr/bin/time ./trou_noir_FP32 $SIZE 1 12 10 $RAIE POSITIVE flux_$SIZE_$RAIE_$HOST.pgm z_$SIZE_$RAIE_$HOST.pgm >>$LOGFILE 2>&1 
done
RAIE=MONOCHROMATIQUE
echo Output stored in $LOGFILE
echo -e "Experience : $RAIE" >>$LOGFILE
echo -e "Experience : $RAIE"
seq $SEQ | while read POWER ; do
    SIZE=$((2**$POWER))
    seq 1 1 10 | xargs -I TOTO  /usr/bin/time ./trou_noir_FP32 $SIZE 1 12 10 $RAIE POSITIVE flux_$SIZE_$RAIE_$HOST.pgm z_$SIZE_$RAIE_$HOST.pgm >>$LOGFILE 2>&1 
done
