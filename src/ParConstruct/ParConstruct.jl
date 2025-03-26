

module ParConstruct


using MPI
using PetscWrap
using FastGaussQuadrature
using FFTW
using JLD2 #perhaps the gather surfs should be somewhere else


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
export gather_surfs


end 





