

#so repl now works on Gadi, and we have confirmed that MID will run as expected.
#TODO
#Make sure we can submit a job with MID, ie check installation location of julia etc! this works!!!
#Make sure MIDParallel will work in login node #this works, now using cpardiso as the main solver as that is intel and pre-installed!
#get MIDParallel to work in job. #Works!!!, data managment and method of doing it is cooked though!

using MIDParallel
#using Plots; plotlyjs() #having this here, and installed in the global environment tricks it into using plotlyjs for interactive plots. This is an awful solution.
#also gives some fkn warning, I think becuase MID doesn't have PlotlyJS.


#simplest case for running jobs with gadi. Bash script will create appropriate inputs, then this will read the inputs and solve.
#not sure if this should be contained within MIDParallel somehow???
dir = ARGS[1]
freq = ARGS[2]

if length(ARGS) > 2
    nev = ARGS[3]
else
    nev = 200
end

#increased number of evals.
#should have that as an input tbh!
par_spectrum_from_file(dir=dir, freq=parse(Float64, freq), nev=nev)
