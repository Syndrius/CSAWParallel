: '
Example script for running CSAWParallel examples in parallel.

This assumes inputs have already been written to file using the Julia REPL.
This file is set up to run from the examples/ subdirectory.
Simply replace these with the appropriate paths.
'

#file created to call par_compute_spectrum
file_to_run=../spectrum_from_file.jl

#where the inputs are written to.
#examples are set up to write into test/data
data_dir=../test/data/

#for qfm cases the surf directory must be specified.
#this points to the example surfaces used for testing
surf_dir=../test/data/benchmark_surfaces.jld2

#Setting the project to CSAWParallel can increase performance
#as global julia environment is not loaded.
PROJ=../

#we can also specify a path to a precompiled .so file 
#precomppilation can improve performance
#SO=/path/to/compiled_file.so

#we can add any extra arguments for slepc via the command line
SLEPC_ARGS=""

#basic execution
mpiexec -n 2 julia --project=$PROJ $file_to_run $data_dir $SLEPC_ARGS
#execution with qfm coordinates
#mpiexec -n 2 julia --project=$PROJ $file_to_run $data_dir $surf_dir $SLEPC_ARGS
#execution with precompilation file
#mpiexec -n 2 julia --project=$PROJ -J$SO $file_to_run $data_dir $SLEPC_ARGS


#we can also post process
#note that this currently requires being run in serial 
julia --project=$PROJ -e "using CSAWParallel; par_post_process(\"$data_dir\")"

