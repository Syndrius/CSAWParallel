"""

"""
module Spectrum

using Printf


using MPI
using PetscWrap
using SlepcWrap


using ChaoticShearAlfvenWaves.WeakForm
using ChaoticShearAlfvenWaves.Io
using ChaoticShearAlfvenWaves.Grids

using ..ParMatrix
using ..Construct
using ..Solve


export par_compute_spectrum


include("QFM.jl") 

export par_compute_spectrum


"""
    par_compute_spectrum(dir::String)

Solves the system PΦ = ω^2QΦ from inputs in the given directory.
"""
function par_compute_spectrum(dir::String)

    uninst_prob, grids, solver = inputs_from_file(dir)

    #this function is quite awkward because jld2 cannot write anonymous functions to file
    #so they have to be recreated.
    prob = WeakForm.inst_problem(uninst_prob.fields, uninst_prob.geo, uninst_prob.flr)

    MPI.Init()

    #if solver isa IntervalSolverT
    #    slepcargs = @sprintf("-eps_interval %f,%f -eps_gen_hermitian -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky st_pc_factor_mat_solver_type superlu_dist -st_mat_superlu_dist_rowperm NOROWPERM", solver.left, solver.right)
    if prob.flr.δ == 0.0 && prob.flr.ρ_i == 0 && prob.flr.δ_e == 0
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
    P, Q = preallocate_matrix(grids)
     
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Constructing...")
    end
    
    #constructs the matrices
    par_construct(P, Q, prob, grids)

    #completes the pestc assembly
    MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY)
    MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY)
    MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY)
    MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        mat_size = matrix_size(grids)
        @printf("Construction of %dx%d matrices complete.\n", mat_size, mat_size)
    end

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Solving...")
    end

    #vectors are used to store the eigenfunctions
    vecr, veci = MatCreateVecs(P)

    #solves the problem
    #and writes the un-processed solutions to file
    nconv = par_solve(P, Q, solver, dir, vecr, veci)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @printf("Solving complete, %d eigenvalues found.\n", nconv)
    end

    #free memory and finalise.
    destroy!(vecr)
    destroy!(veci)
    destroy!(P)
    destroy!(Q)
    SlepcFinalize()
    MPI.Finalize()

end

end
