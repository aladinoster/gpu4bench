#!/bin/sh
gprof -l $1 gmon.out | grep "Pi.c:" | while read LINE; do Line=$(echo $LINE | awk -F: '{ print $2 }' | awk '{ print $1 }') ; Time=$(echo $LINE | awk '{ print $1 }') ; echo $Time:$(sed -n ${Line},${Line}p Pi.c) ; done
