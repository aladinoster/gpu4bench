#!/bin/bash
DEVICE=2
GPU=OpenCL
SEQ="6 1 14"
HOST=$(hostname)
DATE=$(date "+%Y%m%d-%H%M")
LOGFILE=TrouNoir-${HOST}-${DATE}.log
echo >$LOGFILE
[ ! -z "$1" ] && DEVICE="$1"
[ ! -z "$2" ] && GPU="$2"
LINE=BB
METHOD=TrajectoPixel
echo Output stored in $LOGFILE
echo -e "Experience : $LINE $METHOD" >>$LOGFILE
echo -e "Experience : $LINE $METHOD $DEVICE $GPU"
seq $SEQ | while read POWER ; do
    SIZE=$((2**$POWER))
    seq 1 1 10 | xargs -I TOTO /usr/bin/time python ./TrouNoir.py -d $DEVICE -g $GPU -n -b -s $SIZE -o $METHOD >>$LOGFILE 2>&1 
done
METHOD=TrajectoCircle
echo -e "Experience : $LINE $METHOD" >>$LOGFILE
echo -e "Experience : $LINE $METHOD $DEVICE $GPU"
seq $SEQ | while read POWER ; do
    SIZE=$((2**$POWER))
    seq 1 1 10 | xargs -I TOTO /usr/bin/time python ./TrouNoir.py -d $DEVICE -g $GPU -n -b -s $SIZE -o $METHOD >>$LOGFILE 2>&1 
done
METHOD=EachPixel
echo -e "Experience : $LINE $METHOD" >>$LOGFILE
echo -e "Experience : $LINE $METHOD $DEVICE $GPU"
seq $SEQ | while read POWER ; do
    SIZE=$((2**$POWER))
    seq 1 1 10 | xargs -I TOTO /usr/bin/time python ./TrouNoir.py -d $DEVICE -g $GPU -n -b -s $SIZE -o $METHOD >>$LOGFILE 2>&1 
done
LINE=MONO
METHOD=TrajectoPixel
echo -e "Experience : $LINE $METHOD" >>$LOGFILE
echo -e "Experience : $LINE $METHOD $DEVICE $GPU"
seq $SEQ | while read POWER ; do
    SIZE=$((2**$POWER))
    seq 1 1 10 | xargs -I TOTO /usr/bin/time python ./TrouNoir.py -d $DEVICE -g $GPU -n -s $SIZE -o $METHOD >>$LOGFILE 2>&1 
done
METHOD=TrajectoCircle
echo -e "Experience : $LINE $METHOD" >>$LOGFILE
echo -e "Experience : $LINE $METHOD $DEVICE $GPU"
seq $SEQ | while read POWER ; do
    SIZE=$((2**$POWER))
    seq 1 1 10 | xargs -I TOTO /usr/bin/time python ./TrouNoir.py -d $DEVICE -g $GPU -n -s $SIZE -o $METHOD >>$LOGFILE 2>&1 
done
METHOD=EachPixel
echo -e "Experience : $LINE $METHOD" >>$LOGFILE
echo -e "Experience : $LINE $METHOD $DEVICE $GPU"
seq $SEQ | while read POWER ; do
    SIZE=$((2**$POWER))
    seq 1 1 10 | xargs -I TOTO /usr/bin/time python ./TrouNoir.py -d $DEVICE -g $GPU -n -s $SIZE -o $METHOD >>$LOGFILE 2>&1 
done
LINE=BB
METHOD=EachCircle
echo -e "Experience : $LINE $METHOD" >>$LOGFILE
echo -e "Experience : $LINE $METHOD $DEVICE $GPU"
seq $SEQ | while read POWER ; do
    SIZE=$((2**$POWER))
    seq 1 1 10 | xargs -I TOTO /usr/bin/time python ./TrouNoir.py -d $DEVICE -g $GPU -n -s $SIZE -o $METHOD >>$LOGFILE 2>&1 
done
LINE=MONO
METHOD=EachCircle
echo -e "Experience : $LINE $METHOD" >>$LOGFILE
echo -e "Experience : $LINE $METHOD $DEVICE $GPU"
seq $SEQ | while read POWER ; do
    SIZE=$((2**$POWER))
    seq 1 1 10 | xargs -I TOTO /usr/bin/time python ./TrouNoir.py -d $DEVICE -g $GPU -n -s $SIZE -o $METHOD >>$LOGFILE 2>&1 
done
