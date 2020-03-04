#!/bin/bash

variables=( "$@" ) #which variables you want to extract; must be some
                   #combination of "energy", "pressure" and "temperature"

for v in "${variables[@]}"
do
  if [ "$v" == "energy" ]; then
    grep "e_b" REPORT | awk '{ print $2 }' > energy.dat
  fi  
  if [ "$v" == "temperature" ]; then
    grep "t_b" REPORT | awk '{ print $5 }' > temperature.dat
  fi
  if [ "$v" == "pressure" ]; then
    grep "total pressure" OUTCAR | awk '{ print $4 }' > pressure.dat
  fi
done

