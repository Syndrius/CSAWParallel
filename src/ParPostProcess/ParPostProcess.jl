
#directly converting petsc vectors to julia arrays is bugged
#so post processing is done by first writing the petsc solution to hdf5 files
#then later these are read in and processed in the same way as the serial case.
module ParPostProcess

using Printf
using JLD2

using MPI
using PetscWrap
using SlepcWrap
using LinearAlgebra

using MID.PostProcessing
using MID.Io
using MID.Structures



include("PetscToFile.jl")

export par_sols_to_file
export write_evals
export write_errs

include("ProcessFromFile.jl")

export par_post_process


end
