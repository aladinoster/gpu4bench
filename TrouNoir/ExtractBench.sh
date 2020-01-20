#!/bin/bash
egrep '(^Experience|^Size\ |^Elapsed|^Compute)' $1 | tr '\n' '@' | sed -e "s/Experience/\nExperience/g" | sed -e "s/Size/\nSize/g" | awk '{ print $2" "$3" "$4" "$6" "$9 }'  | tr -d ":" | sed -e "s/Compute\ //g" | sed -e "s/Elapsed\ //g" | tr "@" " " | sed -e "s/Time\ //g" | sed -e "s/^\ //g"
