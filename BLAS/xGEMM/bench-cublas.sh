#!/bin/bash

DIR=$(pwd)/OUT

BENCH=lesson23

NUMBER=10

SIZE=30

STEP=1

MAX=3000

FORMAT=SP

OUTSP_THUNKING=$DIR/${BENCH}_${FORMAT}_thunking.out
OUTSP_CUBLAS=$DIR/${BENCH}_${FORMAT}_cublas.out

echo > $OUTSP_THUNKING
echo > $OUTSP_CUBLAS

while [ $SIZE -le $MAX ]
do
        
    FORMAT=SP

    THUNKING=$(./${BENCH}_${FORMAT}_thunking $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
    CUBLAS=$(./${BENCH}_${FORMAT}_cublas $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | tr "\n" " " | awk '{ print  $5"\t"$1"\t"$3 }')
    
    echo -e $SIZE"\t"$THUNKING >> $OUTSP_THUNKING
    echo -e $SIZE"\t"$CUBLAS >> $OUTSP_CUBLAS

    SIZE=$(($SIZE+$STEP))

done
