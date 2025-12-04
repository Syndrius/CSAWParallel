"""
Solves the generalised eigenvalue problem PΦ=ω^2QΦ in parallel using SLEPc.
Solutions are written to file as they are found.
"""
module Solve


using MPI
using PetscWrap
using SlepcWrap


using ChaoticShearAlfvenWaves.Solve


using ..PostProcess


include("SliceSolve.jl")

include("ShiftInvertSolve.jl")

#work in progres...
#incldue("IntervalSolve.jl")

export par_solve



end 
