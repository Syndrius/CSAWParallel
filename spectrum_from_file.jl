"""
File to be called called in parallel that simply calls par_compute_spectrum from command line args.
"""
using CSAWParallel

#first command line arg is the directory where inputs are stored
dir = ARGS[1]

#checks if qfm surfaces have been passed in.
if length(ARGS) > 1 && occursin("jld2", ARGS[2])
    surfs_dir = ARGS[2]

    #removes these two args so that any more slepc args can be read from the command line args starting from one
    deleteat!(ARGS, [1, 2])
    par_compute_spectrum(dir, surfs_dir)
else
    deleteat!(ARGS, 1)
    par_compute_spectrum(dir)
end

