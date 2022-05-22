#!/bin/bash

SIZES=( 16 32 64 128 256)

STEP=( 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 )
TC=2.26919

for L in "${SIZES[@]}"
do  
    TSTART=$( echo "scale=5; ${TC} - 1.0 / $L;" | bc)
    for I in "${STEP[@]}"
    do
        T=$( echo "scale=5; ${TSTART} + 0.2 * $I / $L;" | bc)
        ./ising L=$L nblock=50 T=$T run
    done

done
