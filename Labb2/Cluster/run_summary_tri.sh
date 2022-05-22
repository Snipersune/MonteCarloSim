#!/bin/bash

SIZES=( 016 024 032 048 064 096 128 )

for L in "${SIZES[@]}"
do  
    ls data/tri_${L}* | ./summary > data/sum_tri_${L}.txt
done

