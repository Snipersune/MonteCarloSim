#!/bin/bash

NMAX=200

N=1
while [ $N -lt 101 ]
do  
    ./rwalk N=$N W=500000 run
    N=$(( $N + 1 ))
done

while [ $N -lt $(( ${NMAX} + 1)) ]
do  
    ./rwalk N=$N W=5000000 run
    N=$(( $N + 1 ))
done