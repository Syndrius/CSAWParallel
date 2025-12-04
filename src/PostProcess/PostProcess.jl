"""
Module for post processing the solutions found in parallel.
Notably the current implementation of PetscWrap.jl does not convert PETSc arrays to Julia arrays, so the raw output is written to file and the post-processing must be done after, in serial.
"""
module PostProcess


using Printf
using JLD2
using MPI
using PetscWrap
using SlepcWrap
using LinearAlgebra


using ChaoticShearAlfvenWaves.PostProcessing
using ChaoticShearAlfvenWaves.Io
using ChaoticShearAlfvenWaves.Structures
using ChaoticShearAlfvenWaves.Grids


include("PetscToFile.jl")

export par_sols_to_file
export write_evals
export write_errs


include("ProcessFromFile.jl")

export par_post_process


end
