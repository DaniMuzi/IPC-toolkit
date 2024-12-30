#!/bin/bash

fname=$1

./generator $fname > check.txt 2>&1 &
wait

if [ -s check.txt ]; then
    python3 lattice_generator.py $fname > check2.txt 2>&1 &
    wait

    if [ -s check2.txt ]; then
      echo "The density is too high even for an FCC lattice"
      rm check2.txt
    fi
fi
rm check.txt
