"""
"""
module Solve

using MPI
using PetscWrap
using SlepcWrap

using MID.Solve
using ..PostProcess


include("SliceSolve.jl")

include("ShiftInvertSolve.jl")

#work in progres...
#incldue("IntervalSolve.jl")

export par_solve



end 
