#!/bin/bash

EXE=Pi_MPI_LONG
ITERATIONS=10000000000000
DATE=$(date "+%Y%m%d%H%M")
MyTIME=/usr/bin/time
export TIME="%U %S %e %P %X %D %K %M %I %O %F %R %W %c %w %r %s"

REPEAT=100
PROCESS=326

[ ! $1 == '' ] && EXE=$1
[ ! $2 == '' ] && ITERATIONS=$2
[ ! $3 == '' ] && PROCESS=$3

LOGFILE=${EXE}_${HOSTNAME}_${ITERATIONS}_${DATE}.log

> $LOGFILE
for p in $(seq $PROCESS -1 1)
do
    echo -e "Process $p" >>$LOGFILE
    echo -ne "Start $EXE with $ITERATIONS and $p : "

    for i in $(seq 1 1 $REPEAT)
    do 
	echo -ne "$i "
	$MyTIME mpirun.openmpi -np $p -mca btl self,openib,sm -hostfile /etc/clusters/r410.nodes -loadbalance hwloc-bind -p pu:0-7 ./$EXE $ITERATIONS >>$LOGFILE 2>&1
    done
    echo
done
