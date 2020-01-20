#!/bin/bash

BENCH=xTRSV

NUMBER=1000

SIZE=100

MAX=1000

FORMAT=DP

OUT_CBLAS=/tmp/${BENCH}_${FORMAT}_cblas.out
OUT_FBLAS=/tmp/${BENCH}_${FORMAT}_fblas.out
OUT_GSL=/tmp/${BENCH}_${FORMAT}_gsl.out
OUT_THUNKING=/tmp/${BENCH}_${FORMAT}_thunking.out
OUT_CUBLAS=/tmp/${BENCH}_${FORMAT}_cublas.out
OUT_GOTOBLAS=/tmp/${BENCH}_${FORMAT}_gotoblas.out
OUT_ACML=/tmp/${BENCH}_${FORMAT}_acml.out

echo > $OUT_CBLAS
echo > $OUT_FBLAS
echo > $OUT_GSL
echo > $OUT_THUNKING
echo > $OUT_CUBLAS
echo > $OUT_GOTOBLAS
echo > $OUT_ACML

while [ $SIZE -le $MAX ]
do

    CUBLAS=$(./${BENCH}_${FORMAT}_cublas $SIZE $NUMBER | grep Duration | awk -F: '{ print $2 }' | tr "\n" " " | awk '{ print  $5"\t"$1"\t"$3 }')

    CBLAS=$(./${BENCH}_${FORMAT}_cblas $SIZE $NUMBER | grep Duration | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
    FBLAS=$(./${BENCH}_${FORMAT}_fblas $SIZE $NUMBER | grep Duration | awk -F: '{ print $2 }' | awk '{ print  $1 }')

    GSL=$(./${BENCH}_${FORMAT}_gsl $SIZE $NUMBER | grep Duration | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
    THUNKING=$(./${BENCH}_${FORMAT}_thunking $SIZE $NUMBER | grep Duration | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
    GOTOBLAS=$(./${BENCH}_${FORMAT}_gotoblas $SIZE $NUMBER | grep Duration | awk -F: '{ print $2 }' | tr "\n" " " | awk '{ print  $5"\t"$1"\t"$3 }')
    
    ACML=$(./${BENCH}_${FORMAT}_acml $SIZE $NUMBER | grep Duration | awk -F: '{ print $2 }' | tr "\n" " " | awk '{ print  $5"\t"$1"\t"$3 }')
    
    echo -e $SIZE"\t"$CBLAS >> $OUT_CBLAS
    echo -e $SIZE"\t"$FBLAS >> $OUT_FBLAS
    echo -e $SIZE"\t"$GSL >> $OUT_GSL
    echo -e $SIZE"\t"$THUNKING >> $OUT_THUNKING
    echo -e $SIZE"\t"$CUBLAS >> $OUT_CUBLAS
    echo -e $SIZE"\t"$GOTOBLAS >> $OUT_GOTOBLAS
    echo -e $SIZE"\t"$ACML >> $OUT_ACML

    SIZE=$(($SIZE+100))
    
done
