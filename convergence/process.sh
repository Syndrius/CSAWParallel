#!/usr/bin/env bash
#DIR=/Users/matt/phd/MIDParallel/data/convergence/
DIR=/scratch/y08/mt3516/Helmholtz/test/
source ~/module_files/mid_modules.sh
for d in $DIR*/; do
    echo "$d"
    julia -e 'using MIDParallel; par_post_process("'"${d}"'")'
done
