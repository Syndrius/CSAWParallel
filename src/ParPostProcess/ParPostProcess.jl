"""
Module for post processing the solutions found in parallel.
Notably there is a bug converting the Petsc Vec's into Julia arrays meaning this must be done after each run, reading the results from file.
"""
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
