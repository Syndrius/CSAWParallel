
#stupid name, but matrix is already used...
module ParMatrix

#inconsist with the rest of the module!
using MID
using MID.Structures
using MID.Basis
using MID.Indexing

using MPI
using PetscWrap


include("NonZeroInds.jl")

#this could perhaps be split up.
include("PreAllocate.jl")

export preallocate_matrix

include("MatrixToGrid.jl")

export matrix_to_grid

end
