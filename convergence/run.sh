#!/usr/bin/env bash
#PBS -q normal
#PBS -l walltime=00:30:00
#PBS -l ncpus=17
#PBS -l mem=180Gb
#PBS -l storage=gdata/y08
#PBS -N fff_helm

#DIR=/Users/matt/phd/MIDParallel/data/convergence/
DIR=/scratch/y08/mt3516/Helmholtz/test/
source ~/module_files/mid_modules.sh
for d in $DIR*/; do
    echo "$d"
    mpiexec -n 17 julia -e 'using MIDParallel; par_spectrum_from_file(dir="'"${d}"'")'
done
