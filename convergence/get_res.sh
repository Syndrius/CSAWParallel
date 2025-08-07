#!/usr/bin/env bash
#PBS -q normal
#PBS -l walltime=02:30:00
#PBS -l ncpus=1
#PBS -l mem=30Gb
#PBS -l storage=gdata/y08
#PBS -N res

source ~/module_files/mid_modules.sh
julia /home/149/mt3516/island_damping/MIDParallel/convergence/job_converg.jl
