"""

Module for constructing the P and Q matrices in parallel. Functions are very similar to CSAW, except grid is split across cores and petsc matrices are constructed instead.
"""
module Construct


using MPI
using PetscWrap
using FFTW


using ChaoticShearAlfvenWaves.Structures
using ChaoticShearAlfvenWaves.Geometry
using ChaoticShearAlfvenWaves.Fields
using ChaoticShearAlfvenWaves.Basis
using ChaoticShearAlfvenWaves.WeakForm
using ChaoticShearAlfvenWaves.Grids
using ChaoticShearAlfvenWaves.QFM
using ChaoticShearAlfvenWaves.Integration


using ..ParMatrix


include("FSS.jl")
include("FFS.jl")
include("FFF.jl")

export par_construct


end 





