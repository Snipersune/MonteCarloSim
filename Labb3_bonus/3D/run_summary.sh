#!/bin/bash

SIZES=( 012 016 024 032 048 )

for L in "${SIZES[@]}"
do 
    ls data/${L}* | ./summary > sum/${L}.txt
done