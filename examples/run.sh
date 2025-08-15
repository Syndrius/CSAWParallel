#!/bin/bash
#PBS -q normal
#PBS -l walltime=01:30:00
#PBS -l ncpus=48
#PBS -l mem=180Gb
#PBS -l storage=gdata/y08
#PBS -N superlu_cholesky

#need the last line above so that the job can have access to gdata

source /home/149/mt3516/module_files/mid_modules.sh

#args are not being read, none of this worked
#interval case is args indep but just segfaulted for some reason.
#all very nice.

#default :
#SLEPC_ARGS=-eps_view -eps_gen_hermitian -st_type sinvert -eps_nev 400
#

#superlu_cholesky
#SLEPC_ARGS=-eps_view -eps_gen_hermitian -st_type sinvert -eps_nev 400 -st_pc_type cholesky -st_pc_factor_mat_solver_type superlu_dist
#

#restriction on mpd -> taken to be very small to extremise the example!
#note value of ncv must not be larger than nev+mpd
#SLEPC_ARGS=-eps_view -eps_gen_hermitian -st_type sinvert -eps_nev 400 -st_pc_type cholesky -st_pc_factor_mat_solver_type superlu_dist -eps_mpd 25 -eps_ncv 800

#I node limit :
#SLEPC_ARGS=-eps_view -eps_gen_hermitian -st_type sinvert -eps_nev 400 -mat_inode_limit 500

#Interval :
# -> use the prebuild args, as they are quite specific.

file_to_run=/home/149/mt3516/island_damping/MIDParallel/spectrum_from_file.jl

#data_dir=/scratch/y08/mt3516/slepc_test/superlu_cholesky/
data_dir=/scratch/y08/mt3516/test/

surf_dir=/scratch/y08/mt3516/qfm/surfaces/k05_surfs.jld2

PROJ=/home/149/mt3516/island_damping/MIDParallel


SO=/home/149/mt3516/island_damping/compilation/MIDParallel.so

SLEPC_ARGS="-eps_view -eps_gen_hermitian -st_type sinvert -eps_nev 50"
#SLEPC_ARGS="-eps_view -eps_gen_hermitian -st_type sinvert -eps_nev 50 -st_pc_type cholesky -st_pc_factor_mat_solver_type superlu_dist"
#SLEPC_ARGS="-eps_view -eps_gen_hermitian -st_type sinvert -eps_nev 50 -st_pc_type cholesky -st_pc_factor_mat_solver_type superlu_dist -eps_mpd 25 -eps_ncv 60"
#the i-node limit doesn't do anything!
#SLEPC_ARGS="-eps_view -eps_gen_hermitian -st_type sinvert -eps_nev 400 -mat_inode_limit 500"
#SLEPC_ARGS="-eps_view -eps_harmonic -eps_gen_hermitian -eps_nev 50"

#mpiexec -n 48 julia --project=$PROJ -J$SO $file_to_run $data_dir $surf_dir $SLEPC_ARGS
#mpiexec -n 2 julia --project=$PROJ -J$SO $file_to_run $data_dir $SLEPC_ARGS
mpiexec -n 4 julia --project=$PROJ $file_to_run $data_dir $SLEPC_ARGS
#mpiexec -n 2 julia --project=$PROJ $file_to_run $data_dir

#for interval test!
#mpiexec -n 48 julia --project=$PROJ -J$SO $file_to_run $data_dir $surf_dir

#echo mpiexec -n 48 julia --project=$PROJ -J$SO $file_to_run $data_dir $surf_dir

