
#file for running MIDParallel from file.
#typical use case write the inputs to file then uses this file for execution.

#the existance of this file is pretty fkn annoying.
#unsure how else to implement this tbh.
using MIDParallel

#first command line arg is the directory where inputs are stored
dir = ARGS[1]

#feel like this checking could be done within MIDParallel
#but maybe this is fine.
if length(ARGS) > 1 && occursin("jld2", ARGS[2])
    surfs_dir = ARGS[2]

    #removes these two args so that any more slepc args can be read from the command line args starting from one
    deleteat!(ARGS, [1, 2])
    par_compute_spectrum(dir, surfs_dir)
else
    deleteat!(ARGS, 1)
    par_compute_spectrum(dir)
end

