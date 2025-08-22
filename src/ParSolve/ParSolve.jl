"""
Module for solving the eigenvalue equation, Wϕ = ω^2Iϕ, in parallel.
Uses Slepc to solve the problem, default case uses the Krylov Schur for a specific target frequency.
Slepc arguments can be passed in as command line arguments to modify the process.
"""
module ParSolve

using Printf
using JLD2 


using MPI
using PetscWrap
using SlepcWrap


using MID.Structures
using MID.WeakForm
using MID.Basis
using MID.Io
using MID.QFM
using ..ParMatrix
using ..ParConstruct
using ..ParPostProcess


export par_compute_spectrum
export qfm_compute_spectrum
export par_spectrum_from_file
export qfm_spectrum_from_file


include("ShiftInvertSolve.jl")
include("SliceSolve.jl")
include("IntervalSolve.jl")

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

    if solver isa IntervalSolverT
        slepcargs = @sprintf("-eps_interval %f,%f -eps_gen_hermitian -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky st_pc_factor_mat_solver_type superlu_dist -st_mat_superlu_dist_rowperm NOROWPERM", solver.left, solver.right)
    elseif prob.flr.δ == 0.0 && prob.flr.ρ_i == 0 && prob.flr.δ_e == 0
        #sets the solver to hermitian, adds a few basic default arguments.
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_hermitian", solver.nev)
    else
        #otherwise problem is not Hermitian.
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_non_hermitian", solver.nev)
    end

    #adds any extra slepc args from the command line
    for i in ARGS
        slepcargs *= " " * i
    end

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Initialising slepc with:")
        display(slepcargs)
    end
    #this also inits petsc
    SlepcInitialize(slepcargs)
    
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Preparing Matrices...")
    end
    
    #preallocate the matrix memory.
    W, I = preallocate_matrix(grids)
     
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Constructing...")
    end
    
    #constructs the matrices
    par_construct(W, I, prob, grids)

    #completes the pestc assembly
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

    #solves the problem
    #and writes the un-processed solutions to file
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
    
    if solver isa IntervalSolverT
        slepcargs = @sprintf("-eps_interval %f,%f -eps_gen_hermitian -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky st_pc_factor_mat_solver_type superlu_dist -st_mat_superlu_dist_rowperm NOROWPERM", solver.left, solver.right)
    elseif prob.flr.δ == 0.0 && prob.flr.ρ_i == 0 && prob.flr.δ_e == 0
        #sets the solver to hermitian, adds a few basic default arguments.
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_hermitian", solver.nev)
    else
        #otherwise problem is not Hermitian.
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_non_hermitian", solver.nev)
    end

    #adds any extra slepc args from the command line
    for i in ARGS
        slepcargs *= " " * i
    end

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Initialising slepc with:")
        display(slepcargs)
    end

    #this also inits petsc
    SlepcInitialize(slepcargs)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Preparing Matrices...")
    end
    
    #preallocate the matrix memory.
    W, I = preallocate_matrix(grids)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Constructing...")
    end
    
    #constructs the matrices
    par_construct(W, I, prob, grids, surfs)

    #completes the pestc assembly
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

    #solves the problem
    #and writes the un-processed solutions to file
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

    prob, grids, solver = inputs_from_file(dir=dir)

    par_compute_spectrum(prob=prob, grids=grids, solver=solver, dir=dir)

end

"""
    par_spectrum_from_file(; dir::String, target_freq::Float64, nev=100::Int64)

Computes the spectrum in parallel from inputs stored in files in the given directory.
Case for qfm surfaces as well, subject to lots of change!
"""
function qfm_spectrum_from_file(; dir::String, qfm_surfs::String)

    prob, grids, solver = inputs_from_file(dir=dir)

    surfs = load_object(qfm_surfs)

    qfm_compute_spectrum(prob=prob, grids=grids, solver=solver, surfs=surfs, dir=dir)

end

end
