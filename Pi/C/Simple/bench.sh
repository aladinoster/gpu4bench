#!/bin/bash

EXE=Pi_LONG
ITERATIONS=1000000000
TIME=time

REPEAT=10

[ ! $1 == '' ] && EXE=$1
[ ! $2 == '' ] && ITERATIONS=$2

LOGFILE=${EXE}_${HOSTNAME}_${ITERATIONS}.log

echo -ne "Start $EXE with $ITERATIONS : "
> $LOGFILE
i=1
while [ $i -le $REPEAT ]
do 
    echo -ne "$i "
    $TIME ./$EXE $ITERATIONS >> $LOGFILE 2>&1 
    i=$(($i+1))
done
echo " done."
