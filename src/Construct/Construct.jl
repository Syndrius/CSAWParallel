"""

Module for constructing the W and I matrices in parallel. Functions are very similar to MID, however, petsc matrices are constructed instead.
"""
module Construct


using MPI
using PetscWrap
using FFTW

using MID.Structures
using MID.Geometry
using MID.Fields
using MID.Basis
using MID.WeakForm
using MID.Grids
using MID.QFM
using MID.Integration


using ..ParMatrix


include("FSS.jl")
include("FFS.jl")
include("FFF.jl")

export par_construct


end 





