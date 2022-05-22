#!/bin/bash

SIZES=( 016 032 064 128 256 )


for L in "${SIZES[@]}"
do  
    ls data/${L}* | ./summary > data/sum_${L}.txt
done

