#!/bin/bash


N=1
while [ $N -lt 21 ]
do  
    ./rwalk N=$N W=1000000 run
    N=$(( $N + 1 ))
done