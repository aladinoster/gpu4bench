#!/bin/bash

DIR=$(pwd)/OUT-$(date "+%Y%m%d")-${HOSTNAME}

echo $DIR

[ ! -d $DIR ] && mkdir -p $DIR

BENCH=xGEMM

NUMBER=2

STEP=1000

START=1000

ENDS=16000
ENDD=16000

FORMAT=SP

OUTSP_CBLAS=$DIR/${BENCH}_${FORMAT}_cblas.out
OUTSP_FBLAS=$DIR/${BENCH}_${FORMAT}_fblas.out
OUTSP_GSL=$DIR/${BENCH}_${FORMAT}_gsl.out
OUTSP_THUNKING=$DIR/${BENCH}_${FORMAT}_thunking.out
OUTSP_CUBLAS=$DIR/${BENCH}_${FORMAT}_cublas.out
OUTSP_OPENBLAS=$DIR/${BENCH}_${FORMAT}_openblas.out
OUTSP_ACML=$DIR/${BENCH}_${FORMAT}_acml.out

FORMAT=DP

OUTDP_CBLAS=$DIR/${BENCH}_${FORMAT}_cblas.out
OUTDP_FBLAS=$DIR/${BENCH}_${FORMAT}_fblas.out
OUTDP_GSL=$DIR/${BENCH}_${FORMAT}_gsl.out
OUTDP_THUNKING=$DIR/${BENCH}_${FORMAT}_thunking.out
OUTDP_CUBLAS=$DIR/${BENCH}_${FORMAT}_cublas.out
OUTDP_OPENBLAS=$DIR/${BENCH}_${FORMAT}_openblas.out
OUTDP_ACML=$DIR/${BENCH}_${FORMAT}_acml.out

echo > $OUTSP_CBLAS
echo > $OUTSP_FBLAS
echo > $OUTSP_GSL
echo > $OUTSP_THUNKING
echo > $OUTSP_CUBLAS
echo > $OUTSP_ACML
echo > $OUTSP_OPENBLAS

echo > $OUTDP_CBLAS
echo > $OUTDP_FBLAS
echo > $OUTDP_GSL
echo > $OUTDP_THUNKING
echo > $OUTDP_CUBLAS
echo > $OUTDP_ACML
echo > $OUTDP_OPENBLAS

SIZE=$START

while [ $SIZE -le $ENDS ]
do
        
    FORMAT=SP

    THUNKING=$(./${BENCH}_${FORMAT}_thunking $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
    CUBLAS=$(./${BENCH}_${FORMAT}_cublas $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | tr "\n" " " | awk '{ print  $5"\t"$1"\t"$3 }')
    
#    CBLAS=$(./${BENCH}_${FORMAT}_cblas $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
#    FBLAS=$(./${BENCH}_${FORMAT}_fblas $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')

#    GSL=$(./${BENCH}_${FORMAT}_gsl $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
#    OPENBLAS=$(./${BENCH}_${FORMAT}_openblas $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
#    ACML=$(./${BENCH}_${FORMAT}_acml $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
    echo -e $SIZE"\t"$THUNKING >> $OUTSP_THUNKING
    echo -e $SIZE"\t"$CUBLAS >> $OUTSP_CUBLAS
    echo -e $SIZE"\t"$CBLAS >> $OUTSP_CBLAS
    echo -e $SIZE"\t"$FBLAS >> $OUTSP_FBLAS
    echo -e $SIZE"\t"$GSL >> $OUTSP_GSL
    echo -e $SIZE"\t"$OPENBLAS >> $OUTSP_OPENBLAS
    echo -e $SIZE"\t"$ACML >> $OUTSP_ACML

    SIZE=$(($SIZE+$STEP))
    
done

SIZE=$START

while [ $SIZE -le $ENDD ]
do

    FORMAT=DP

    THUNKING=$(./${BENCH}_${FORMAT}_thunking $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
    CUBLAS=$(./${BENCH}_${FORMAT}_cublas $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | tr "\n" " " | awk '{ print  $5"\t"$1"\t"$3 }')
    
#    CBLAS=$(./${BENCH}_${FORMAT}_cblas $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')

#    FBLAS=$(./${BENCH}_${FORMAT}_fblas $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')

#    GSL=$(./${BENCH}_${FORMAT}_gsl $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
#    OPENBLAS=$(./${BENCH}_${FORMAT}_openblas $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
#    ACML=$(./${BENCH}_${FORMAT}_acml $SIZE $NUMBER | grep GFlops | awk -F: '{ print $2 }' | awk '{ print  $1 }')
    
    echo -e $SIZE"\t"$CBLAS >> $OUTDP_CBLAS
    echo -e $SIZE"\t"$FBLAS >> $OUTDP_FBLAS
    echo -e $SIZE"\t"$GSL >> $OUTDP_GSL
    echo -e $SIZE"\t"$THUNKING >> $OUTDP_THUNKING
    echo -e $SIZE"\t"$CUBLAS >> $OUTDP_CUBLAS
    echo -e $SIZE"\t"$OPENBLAS >> $OUTDP_OPENBLAS
    echo -e $SIZE"\t"$ACML >> $OUTDP_ACML

    SIZE=$(($SIZE+$STEP))
    
done
