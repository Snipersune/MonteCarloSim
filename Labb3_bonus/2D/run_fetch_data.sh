#!/bin/bash

SIZES=( 32 64 128 256)
PS=( 0.5 0.52 0.54 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64 0.66 0.68 0.7 )

for L in "${SIZES[@]}"
do 
    for P in "${PS[@]}"
    do
        ./perc L=$L N=1000 p=$P run
    done
done