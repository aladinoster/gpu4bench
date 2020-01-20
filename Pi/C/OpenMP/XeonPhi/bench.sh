#!/bin/bash

EXE=Pi_OpenMP_FP32_MWC
#EXE=Pi_OpenMP_FP64_MWC
ITERATIONS=100000000000
MYTIME=/usr/bin/time

REPEAT=10
START=960
END=1920

[ ! $1 == '' ] && EXE=$1
[ ! $2 == '' ] && ITERATIONS=$2
[ ! $3 == '' ] && PROCESS=$3

LOGFILE=${EXE}_${HOSTNAME}_${ITERATIONS}_$(date "+%Y%m%d").log

> $LOGFILE
for p in $(seq $START $END)
do
    export OMP_NUM_THREADS=$p
    echo -e "Process $p" >> $LOGFILE
    echo -ne "Start $EXE with $ITERATIONS and $p : "
    for i in $(seq $REPEAT)
    do 
        echo -ne "$i "
        #$MYTIME hwloc-bind -p pu:1 ./$EXE $ITERATIONS $p >> $LOGFILE 2>&1 
        $MYTIME ./$EXE $ITERATIONS $p >> $LOGFILE 2>&1 
	sleep 10
    done
    echo
done
