
#TODO
#needs to be more general in the final state

file_to_run=/Users/matt/phd/MIDParallel/spectrum_from_file.jl

data_dir=/Users/matt/phd/MIDParallel/data/example/
#data_dir=/Users/matt/phd/MIDParallel/data/mapping/

surf_dir=/Users/matt/phd/MID/test/data/benchmark_surfaces.jld2
#surf_dir=/Users/matt/phd/MIDParallel/data/example/w1_surfaces.jld2

PROJ=/Users/matt/phd/MIDParallel

#if compilation is a go.
#SO=/home/149/mt3516/island_damping/compilation/MIDParallel.so

SLEPC_ARGS=""

#mpiexec -n 2 julia --project=$PROJ -J$SO $file_to_run $data_dir $SLEPC_ARGS
#mpiexec -n 2 julia --project=$PROJ $file_to_run $data_dir $SLEPC_ARGS
mpiexec -n 2 julia --project=$PROJ $file_to_run $data_dir $surf_dir $SLEPC_ARGS
#mpiexec -n 2 julia $file_to_run $data_dir $SLEPC_ARGS


