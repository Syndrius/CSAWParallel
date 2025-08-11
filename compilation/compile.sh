#!/bin/bash
#ideally this file will actually create the folder where data is stored etc
#and create inputs for a qfm test case as well!
#so this doesn't look like it will work on mac due to MPI versions
#hopefully this will be ok on Gadi
#as there we have more clarity on what MPI version is being used everywhere.
DIR=/Users/matt/phd/MIDParallel/compilation/
PROJ=/Users/matt/phd/MIDParallel
mkdir $DIR"data/"

#julia --project=/Users/matt/phd/MIDParallel compile.jl
mpiexec -n 1 julia --project=$PROJ warmup.jl $DIR"data/"

rm -r $DIR"data/"
