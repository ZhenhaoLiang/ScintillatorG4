#!/bin/bash

MC_HOME='.'
for i in $(seq 1 50)
  do
    export Filename='out/'$i
    export Logfile='out/log'$i'.txt'
    $MC_HOME/build/toyMC run.mac $Filename i >$Logfile &
    echo "$i" 
  done

