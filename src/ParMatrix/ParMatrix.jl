"""

Module containing functions to handle the matrices in parallel. This includes preallocation and indexing across cores.
"""
module ParMatrix


using ChaoticShearAlfvenWaves.Structures
using ChaoticShearAlfvenWaves.Grids


using MPI
using PetscWrap


include("NonZeroInds.jl")


include("PreAllocate.jl")

export preallocate_matrix


include("MatrixToGrid.jl")

export matrix_to_grid


end
