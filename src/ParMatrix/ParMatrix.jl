"""

Module containing functions to handle the matrices in parallle. This includes preallocation and indexing across cores.
"""
module ParMatrix


using MID.Structures
using MID.Basis
using MID.Indexing


using MPI
using PetscWrap



include("NonZeroInds.jl")


include("PreAllocate.jl")

export preallocate_matrix



include("MatrixToGrid.jl")

export matrix_to_grid

end
