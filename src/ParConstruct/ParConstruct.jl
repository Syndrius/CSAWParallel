"""

Module for constructing the W and I matrices in parallel. Functions are very similar to MID, however, petsc matrices are constructed instead.
"""
module ParConstruct


using MPI
using PetscWrap
using FastGaussQuadrature
using FFTW

using MID.Structures
using MID.Geometry
using MID.Equilibrium
using MID.Basis
using MID.WeakForm
using MID.Indexing
using MID.QFM
using MID.Integration


using ..ParMatrix


include("FSS.jl")
include("FFS.jl")
include("FFF.jl")

export par_construct


end 





