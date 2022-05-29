#!/bin/bash

SIZES=( 032 064 128 256)

for L in "${SIZES[@]}"
do 
    ls data/${L}* | ./summary > sum/${L}.txt
done