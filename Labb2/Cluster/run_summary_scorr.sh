#!/bin/bash

SIZES=( 256 1024 2048 4096 )


for L in "${SIZES[@]}"
do  
    ls data/scorr_${L}* | ./summary > data/sum_scorr_${L}.txt
done

