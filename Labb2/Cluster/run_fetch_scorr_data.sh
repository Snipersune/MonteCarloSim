#!/bin/bash

TC=2.26919

TEMPS=( 2.16 2.20 2.22 2.24 2.26 2.26919 2.28 2.30 2.32 2.36 ) 

for T in "${TEMPS[@]}"
do  
    ./ising L=256 nblock=50 T=$T run
done

./ising L=1024 nblock=20 T=${TC} run
./ising L=2048 nblock=20 T=${TC} run
./ising L=4096 nblock=20 T=${TC} run