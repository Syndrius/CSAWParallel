#!/usr/bin/env bash
DIR=/Users/matt/phd/MIDParallel/data/convergence/
for d in $DIR*/; do
    echo "$d"
    mpiexec -n 2 julia -e 'using MIDParallel; par_spectrum_from_file(dir="'"${d}"'")'
done
