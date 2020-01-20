#!/bin/bash
egrep '(^Experience|Dimension\ |^Elapsed)' $1 | tr '\n' ' ' | sed -e "s/Experience/\nExperience/g" | tr -d "#'" | sed -e "s/Dimension\ de\ limage\ :\ /\n/g" | awk '{ print $1" "$NF }' 
