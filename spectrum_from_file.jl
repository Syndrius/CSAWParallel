

#so repl now works on Gadi, and we have confirmed that MID will run as expected.
#TODO
#Make sure we can submit a job with MID, ie check installation location of julia etc! this works!!!
#Make sure MIDParallel will work in login node #this works, now using cpardiso as the main solver as that is intel and pre-installed!
#get MIDParallel to work in job. #Works!!!, data managment and method of doing it is cooked though!

using MIDParallel
#using Plots; plotlyjs() #having this here, and installed in the global environment tricks it into using plotlyjs for interactive plots. This is an awful solution.
#also gives some fkn warning, I think becuase MID doesn't have PlotlyJS.

#this function should be able to do both normal case and qfm case now.
#simplest case for running jobs with gadi. Bash script will create appropriate inputs, then this will read the inputs and solve.
#not sure if this should be contained within MIDParallel somehow???
dir = ARGS[1]
#freq = ARGS[2]

#second arg is a surface list!
if length(ARGS) > 1 && occursin("jld2", ARGS[2])
    surfs_dir = ARGS[2]

    #removes these two args so that any more slepc args can be read from the command line args starting from one
    deleteat!(ARGS, [1, 2])
    qfm_spectrum_from_file(dir=dir, qfm_surfs=surfs_dir)
else
    deleteat!(ARGS, 1)
    par_spectrum_from_file(dir=dir)
end

#increased number of evals.
#should have that as an input tbh!
#par_spectrum_from_file(dir=dir, target_freq=parse(Float64, freq), nev=nev)
