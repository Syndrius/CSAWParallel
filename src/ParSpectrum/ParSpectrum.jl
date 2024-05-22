
module ParSpectrum

#so this works, obvs without registering MID, we have to add this locally, ie Pkg.add ~/phd/MID or whatever
using MID #must be added locally as we do not have MID on the registary yet.
using MPI
using FFTW
using FastGaussQuadrature
using SparseArrays
using Printf


#note to get these files to work on mac, had to modify the load.jl files in both cases
#mac uses dylib files, while linux uses .so files so PetscWrap and SlepcWrap were unable to find the files.
using PetscWrap
using SlepcWrap


export par_construct_and_solve


include("ParConstruct.jl")

export par_construct


include("ParSolve.jl")

export par_solve


"""
Computes the spectrum in parallel, and writes solution to file.

# Args
- prob::ProblemT Struct containing the functions and parameters that define the problem we are solving
- grids::GridT Grids to solve over.
- σ::Float64 Frequency for 'shift and invert', typically this is the TAE frequency.
- nev::Int64 Number of eigenvalues to solve for.
- dir::String Directory the results are written to.
"""
function par_construct_and_solve(; prob::MID.ProblemT, grids::MID.GridsT, σ=0.0::Float64, nev=20::Int64, dir::String)

    MPI.Init()

    #construct the matrices, each proc has a subset of the coo data.        
    rows, cols, Wdata, Idata = par_construct(prob=prob, grids=grids)

    #compute the matrix size for creating the Petsc matrix.
    n = 2 * grids.rd.N * grids.pmd.count * grids.tmd.count

    #solve the matrix, this writes to file
    par_solve(rows, cols, Wdata, Idata, σ=σ, nev=nev, n=n, dir=dir)

    MPI.Finalize()

end


end