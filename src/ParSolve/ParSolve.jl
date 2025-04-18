"""

Module for solving the eigenvalue equation, Wϕ = ω^2Iϕ, in parallel.
"""
module ParSolve

using Printf
using JLD2 


#note to get these files to work on mac, had to modify the load.jl files in both cases
#mac uses dylib files, while linux uses .so files so PetscWrap and SlepcWrap were unable to find the files.
using MPI
using PetscWrap
using SlepcWrap


using MID.PostProcessing
using MID.Structures
using MID.WeakForm
using MID.Basis
using MID.Io
using MID.QFM
using ..ParMatrix
using ..ParConstruct


export par_compute_spectrum
export qfm_compute_spectrum
export par_spectrum_from_file
export qfm_spectrum_from_file



include("PostProcess.jl")

export par_post_process #ideally this would not be needed but some weird shit is happening.


include("ShiftInvertSolve.jl")


include("SliceSolve.jl")

export par_solve


"""
    par_compute_spectrum(; prob::ProblemT, grids::GridsT, solver::SolverT, dir::String)

Computes the spectrum in parallel, and writes solution to file.

# Args
- prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving.
- grids::GridT - Grids to solve over.
- solver::SolverT - solver struct specifiying how to solve the system.
- dir::String Directory the results are written to.
"""
function par_compute_spectrum(; prob::ProblemT, grids::GridsT, solver::SolverT, dir::String)

    MPI.Init()

    #ideally this should be in the other solver types, even if it is unused in serial, it may be important in slepc.
    if prob.flr.δ == 0.0 && prob.flr.ρ_i == 0 && prob.flr.δ_e == 0
        #sets the solver to hermitian, unsure if it actually matters.
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_hermitian -eps_view", solver.nev) #* evals_str #* efuncs_str 
    else
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_non_hermitian -eps_view", solver.nev) #* evals_str #* efuncs_str 
    end

    #this also inits petsc
    SlepcInitialize(slepcargs)
    
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Preparing Matrices...")
    end
    
    #preallocate the matrix memory.
    #this has an large effect on memory usage.
    W, I = preallocate_matrix(grids)
     
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Constructing...")
    end
    
    par_construct(W, I, prob, grids)
    
    MatAssemblyBegin(W, MAT_FINAL_ASSEMBLY)
    MatAssemblyBegin(I, MAT_FINAL_ASSEMBLY)
    MatAssemblyEnd(W, MAT_FINAL_ASSEMBLY)
    MatAssemblyEnd(I, MAT_FINAL_ASSEMBLY)


    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        mat_size = matrix_size(grids)
        @printf("Construction of %dx%d matrices complete.\n", mat_size, mat_size)
    end

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Solving...")
    end

    #vectors are used to store the eigenfunctions
    #created here for efficiency
    vecr, veci = MatCreateVecs(W)

    nconv = par_solve(W, I, solver, dir, vecr, veci)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @printf("Solving complete, %d eigenvalues found.\n", nconv)
    end

    #free memory and finalise.
    destroy!(vecr)
    destroy!(veci)
    destroy!(W)
    destroy!(I)
    SlepcFinalize()
    MPI.Finalize()

end


"""
    qfm_compute_spectrum(; prob::ProblemT, grids::GridsT, solver::SolverT, dir::String, surfs::Array{QFMSurfaceT})

Constructs the matrices and solves for the case with qfm surfaces.
"""
function qfm_compute_spectrum(; prob::ProblemT, grids::GridsT, solver::SolverT, dir::String, surfs::Array{QFMSurfaceT})

    MPI.Init()
    
    #ideally this should be in the other solver types, even if it is unused in serial, it may be important in slepc.
    if prob.flr.δ == 0.0 && prob.flr.ρ_i == 0 && prob.flr.δ_e == 0
        #sets the solver to hermitian, unsure if it actually matters.
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_hermitian -eps_view", solver.nev) #* evals_str #* efuncs_str 
    else
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_non_hermitian -eps_view", solver.nev) #* evals_str #* efuncs_str 
    end

    #this also inits petsc
    SlepcInitialize(slepcargs)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Preparing Matrices...")
    end
    
    #preallocate the matrix memory.
    #this has an large effect on memory usage.
    W, I = preallocate_matrix(grids)

     
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Constructing...")
    end
    
    par_construct(W, I, prob, grids, surfs)

    
    MatAssemblyBegin(W, MAT_FINAL_ASSEMBLY)
    MatAssemblyBegin(I, MAT_FINAL_ASSEMBLY)
    MatAssemblyEnd(W, MAT_FINAL_ASSEMBLY)
    MatAssemblyEnd(I, MAT_FINAL_ASSEMBLY)


    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        mat_size = matrix_size(grids)
        @printf("Construction of %dx%d matrices complete.\n", mat_size, mat_size)
    end


    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Solving...")
    end

    vecr, veci = MatCreateVecs(W)

    nconv = par_solve(W, I, solver, dir, vecr, veci)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @printf("Solving complete, %d eigenvalues found.\n", nconv)
    end

    #free memory and finalise.
    destroy!(vecr)
    destroy!(veci)
    destroy!(W)
    destroy!(I)
    SlepcFinalize()
    MPI.Finalize()

end


"""
    par_spectrum_from_file(; dir::String, target_freq::Float64, nev=100::Int64)

Computes the spectrum in parallel from inputs stored in files in the given directory.
"""
function par_spectrum_from_file(; dir::String)

    #should only root be doing this?
    prob, grids, solver = inputs_from_file(dir=dir)

    par_compute_spectrum(prob=prob, grids=grids, solver=solver, dir=dir)

end

"""
    par_spectrum_from_file(; dir::String, target_freq::Float64, nev=100::Int64)

Computes the spectrum in parallel from inputs stored in files in the given directory.
Case for qfm surfaces as well, subject to lots of change!
"""
function qfm_spectrum_from_file(; dir::String, qfm_surfs::String)

    #should only root be doing this?
    prob, grids, solver = inputs_from_file(dir=dir)

    surfs = load_object(qfm_surfs)

    qfm_compute_spectrum(prob=prob, grids=grids, solver=solver, surfs=surfs, dir=dir)

end

end
