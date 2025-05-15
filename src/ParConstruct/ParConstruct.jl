"""

Module for constructing the W and I matrices in parallel. Functions are very similar to MID, however, petsc matrices are constructed instead.
Note that this function is still in a bit of a state. Needs lots of comments and organising, but should be only surface level changes.
"""
module ParConstruct


using MPI
using PetscWrap
using FastGaussQuadrature
using FFTW
using JLD2 #perhaps the gather surfs should be somewhere else
using Printf #temp for testing


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



include("QFMSurfaces.jl")

export par_construct_surfaces
export gather_surfs #ideally this will be fixed if we ever require significantly more surfaces


end 





