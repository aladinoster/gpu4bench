#!/bin/bash

EXE=Pi_Pthreads_LONG
ITERATIONS=10000000000
TIME=time

REPEAT=10
PROCESS=16

[ ! $1 == '' ] && EXE=$1
[ ! $2 == '' ] && ITERATIONS=$2
[ ! $3 == '' ] && PROCESS=$3

LOGFILE=${EXE}_${HOSTNAME}_${ITERATIONS}.log

> $LOGFILE
p=1
while [ $p -le $PROCESS ]
do
    echo -e "Process $p" >> $LOGFILE
    echo -ne "Start $EXE with $ITERATIONS and $p : "
    i=1
    while [ $i -le $REPEAT ]
    do 
        echo -ne "$i "
        $TIME ./$EXE $ITERATIONS $p >> $LOGFILE 2>&1 
        i=$(($i+1))
    done
    echo
    p=$(($p+1))
done
