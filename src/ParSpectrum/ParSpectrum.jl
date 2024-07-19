
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
#new version that uses petsc matrices entirely to try and reduce memory usage.
function par_construct_and_solve(; prob::MID.ProblemT, grids::MID.GridsT, σ=0.0::Float64, nev=50::Int64, dir::String)

    target_freq = σ^2 / prob.geo.R0^2 
    MPI.Init()
    #eigenvalues are written in matlab format as this offers greater precision than default, 
    #and is easier to read from file
    evals_str = " -eps_view_values :" * dir * "vals.dat:ascii_matlab"
    #eigenfunctions are written as symmodu, which ignores which part of the eigenfunctions is stored by each proc.
    efuncs_str = " -eps_view_vectors :" * dir * "funcs.dat:ascii_symmodu"

    #these are combined with nev and σ, specifiying the setup for slepc.
    #eps_harmonics means we are searching inside the spectrum.
    #-st_type sinvert
    #slepcargs = @sprintf("-eps_nev %d -eps_target %s -eps_harmonic", nev, σ) * evals_str * efuncs_str
    #-eps_view for solver stuff
    #-memory_view for mem
    #log_view for heaps of petsc info.
    slepcargs = @sprintf("-eps_nev %d -eps_target %s -st_type sinvert -memory_view", nev, target_freq) * evals_str * efuncs_str
    
    #this should also init petsc
    SlepcInitialize(slepcargs)
    #compute the matrix size for creating the Petsc matrix.
    n = matrix_dim(grids) #should work for both grids now!
    #n = 2 * grids.r.N * grids.pmd.count * grids.tmd.count
    W = create_matrix(n, n, auto_setup=true)
    I = create_matrix(n, n, auto_setup=true)

    #construct the matrices, each proc has a subset of the coo data.        
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Constructing...")
    end
    #par_construct(W, I, prob=prob, grids=grids)


    

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Solving...")
    end
    #memory doubling does not seem to be in solve...
    #solve the matrix, this writes to file
    #par_solve(W, I)

    destroy!(W)
    destroy!(I)
    SlepcFinalize()
    MPI.Finalize()

end



"""
Computes the spectrum in parallel, and writes solution to file.

# Args
- prob::ProblemT Struct containing the functions and parameters that define the problem we are solving
- grids::GridT Grids to solve over.
- σ::Float64 Frequency for 'shift and invert', typically this is the TAE frequency.
- nev::Int64 Number of eigenvalues to solve for.
- dir::String Directory the results are written to.
"""
function old_par_construct_and_solve(; prob::MID.ProblemT, grids::MID.GridsT, σ=0.0::Float64, nev=5::Int64, dir::String)

    MPI.Init()

    #construct the matrices, each proc has a subset of the coo data.        
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Constructing...")
    end
    rows, cols, Wdata, Idata = old_par_construct(prob=prob, grids=grids)

    #compute the matrix size for creating the Petsc matrix.
    n = 2 * grids.rd.N * grids.pmd.count * grids.tmd.count

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Solving...")
    end
    #solve the matrix, this writes to file
    old_par_solve(rows, cols, Wdata, Idata, σ=σ, nev=nev, n=n, dir=dir)

    MPI.Finalize()

end


end