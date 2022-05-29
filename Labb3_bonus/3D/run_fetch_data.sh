#!/bin/bash

SIZES=( 12 16 24 32 48 64 96 128 192 )
PS=( 0.29 0.295 0.30 0.305 0.3075 0.31 0.311 0.312 0.313 0.314 0.315 0.3175 0.32 0.325 0.33 )

for L in "${SIZES[@]}"
do 
    for P in "${PS[@]}"
    do
        ./perc L=$L N=10000 p=$P run
    done
done