
#file for running MIDParallel from file.
#typical use case write the inputs to file then uses this file for execution.

using MIDParallel

#first command line arg is the directory where inputs are stored
dir = ARGS[1]

#checks if the second arg is a list of qfm surfaces
#so that qfm case can be run.
if length(ARGS) > 1 && occursin("jld2", ARGS[2])
    surfs_dir = ARGS[2]

    #removes these two args so that any more slepc args can be read from the command line args starting from one
    deleteat!(ARGS, [1, 2])
    qfm_spectrum_from_file(dir=dir, qfm_surfs=surfs_dir)
else
    deleteat!(ARGS, 1)
    par_spectrum_from_file(dir=dir)
end

