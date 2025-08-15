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


#using MID.PostProcessing
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



#include("PostProcess.jl")

#export par_post_process #ideally this would not be needed but some weird shit is happening.


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

    #ideally this should be in the other solver types, even if it is unused in serial, it may be important in slepc.
    #TODO
    #fix the solver to already have this info
    #add option for the diagonal shift
    #that seems to fix the solid vertical lines problem
    #but is a bit arbitrary
    #choosing too large changes the spectrum
    #too small and nothing happens.
    #also as a command line arg it doesn't want 1e-6, that gets read as 1.
    #which may cause problemos for putting in a float tbh.
    #the shift doesn't seem to do anything when it is actually run in parallel
    #yip fkn ee
    #guess we will ignore for now.
    #MUMPS is a bit better, but we are still facing issues
    #we probably need to shift the zero_pivot tolerance if possible!

    #looks like serious memory problems are created with st_pc_type lu.
    #perhaps it defaults to this for low memory cases and defaults to an iterative solver otherwise?
    #looks like we should be doing cholesky factorisation instead of lu, it should be the same but for Hermitian matrices, making it more efficeint.
    #may also want to look at spectrum slicing more seriously. -> this will require a new solver type!
    #think we need to do some testing on a smallish case.

    #SLEPCARGS of interest
    #-eps_true_residual -> compute the residual of the original problem not the shift inverted problem, may be moreaccurate.
    #-eps_mpd may help with memory scale with large number of processes, it looks like there is always a serial step, this restricts the size of that step, which may increase number of iterations but reduce memory
    #-eps_harmonic may be an alterantive to the shift and invert transformation for finding interior evals.
    if solver isa IntervalSolverT
        #this just doesn't work on Gadi, gives seg fault, unsure why, probably not worth trying to fix.
        slepcargs = @sprintf("-eps_interval %f,%f -eps_gen_hermitian -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky st_pc_factor_mat_solver_type superlu_dist -st_mat_superlu_dist_rowperm NOROWPERM", solver.left, solver.right)
    elseif prob.flr.δ == 0.0 && prob.flr.ρ_i == 0 && prob.flr.δ_e == 0
        #sets the solver to hermitian, unsure if it actually matters.
        #slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_hermitian -eps_view -st_pc_type lu -st_pc_factor_mat_solver_type superlu_dist -st_pc_factor_shift_type nonzero -st_pc_factor_shift_amount 0.00001", solver.nev) #* evals_str #* efuncs_str 
        #slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_hermitian -eps_view -st_pc_type lu -st_pc_factor_mat_solver_type superlu_dist -mat_superlu_dist_replacetinypivot -st_pc_factor_shift_type positive_definite", solver.nev) #* evals_str #* efuncs_str 
        #slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_hermitian -eps_view -mat_type mpisbaij", solver.nev) #* evals_str #* efuncs_str 
        #writing the same command twice just overwrites the first
        #I think that works perf for us
        #as we can just have a few of the default ones.
        #but then they can be overwritten externally.
        #slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_hermitian -eps_view -eps_error_relative ::ascii_info_detail -eps_nev 1", solver.nev) #* evals_str #* efuncs_str 
        #default to superlu_dist as that seems to be better for our cases, however this can be overwritten with command line args
        #superlu seems to be extremely more demanding on memory! unsure why, we will stick to default
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_hermitian -st_pc_type lu -st_pc_factor_mat_solver_type superlu_dist", solver.nev)
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_hermitian", solver.nev)
        #slepcargs = @sprintf("-eps_nev %d -eps_gen_hermitian", solver.nev)
        #slepcargs = ""
    else
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_non_hermitian -eps_view", solver.nev) #* evals_str #* efuncs_str 
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_non_hermitian", solver.nev)
        #default to superlu_dist as that seems to be better for our cases, however this can be overwritten with command line args
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_non_hermitian -st_pc_type lu -st_pc_factor_mat_solver_type superlu_dist", solver.nev)
    end

    #display(ARGS)
    for i in ARGS
        #display(i)
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

    #display(MatIsHermitian(I))
    #display(MatIsHermitian(W))

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        mat_size = matrix_size(grids)
        @printf("Construction of %dx%d matrices complete.\n", mat_size, mat_size)
    end

    #=
    viewerW = PetscViewerASCIIOpen(MPI.COMM_WORLD, "/Users/matt/phd/W.dat")
    viewerI = PetscViewerASCIIOpen(MPI.COMM_WORLD, "/Users/matt/phd/I.dat")

    #this doesn't even write the complex part... wot the fek
    MatView(W, viewerW)
    MatView(I, viewerI)

    PetscViewerDestroy(viewerW)
    PetscViewerDestroy(viewerI)
    =#
    

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Solving...")
    end

    #vectors are used to store the eigenfunctions
    #created here for efficiency
    vecr, veci = MatCreateVecs(W)

    nconv = par_solve(W, I, solver, dir, vecr, veci)
    #nconv = 0

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

function MatIsHermitian(mat::PetscMat)

    #so fss is 
    tol = 1e-15 #so similar to other case, Hermitian approximatly, at least for fss.
    #bigger issue may arise when we have a huge number of cores,
    #or when there are large island / or qfm etc
    result = Ref{PetscWrap.PetscBool}(PetscWrap.PETSC_FALSE)

    #this just seems to not work properly, not ideal
    #this is a completly useless check
    #either doesn't work
    #or if we explicity set the matrix to hermitian it just returns true, even if matrix is very non hermitian.
    error = ccall((:MatIsHermitian, PetscWrap.libpetsc), PetscErrorCode, (Ptr{Cvoid}, PetscReal, Ref{PetscWrap.PetscBool}), mat, tol, result)
    #error = ccall((:MatIsSymmetric, PetscWrap.libpetsc), PetscErrorCode, (Ptr{Cvoid}, PetscReal, Ref{PetscWrap.PetscBool}), mat, tol, result)

    #lets try actually getting the transpose and seeing the difference.
    @assert iszero(error)

    return result[]
end


"""
    qfm_compute_spectrum(; prob::ProblemT, grids::GridsT, solver::SolverT, dir::String, surfs::Array{QFMSurfaceT})

Constructs the matrices and solves for the case with qfm surfaces.
"""
function qfm_compute_spectrum(; prob::ProblemT, grids::GridsT, solver::SolverT, dir::String, surfs::Array{QFMSurfaceT})

    MPI.Init()
    
    #ideally this should be in the other solver types, even if it is unused in serial, it may be important in slepc.
    if solver isa IntervalSolverT
        slepcargs = @sprintf("-eps_interval %f,%f -eps_gen_hermitian -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky st_pc_factor_mat_solver_type superlu_dist -st_mat_superlu_dist_rowperm NOROWPERM", solver.left, solver.right)
    elseif prob.flr.δ == 0.0 && prob.flr.ρ_i == 0 && prob.flr.δ_e == 0
        #sets the solver to hermitian, unsure if it actually matters.
        #slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_hermitian -eps_view", solver.nev) #* evals_str #* efuncs_str 
        #slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view -eps_gen_hermitian -eps_view", solver.nev) #* evals_str #* efuncs_str 
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_hermitian", solver.nev)
        #default to superlu_dist as that seems to be better for our cases, however this can be overwritten with command line args
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_hermitian -st_pc_type lu -st_pc_factor_mat_solver_type superlu_dist", solver.nev)
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_hermitian", solver.nev)
    else
        #slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_non_hermitian -eps_view", solver.nev) #* evals_str #* efuncs_str 
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_non_hermitian", solver.nev)
        #default to superlu_dist as that seems to be better for our cases, however this can be overwritten with command line args
        slepcargs = @sprintf("-eps_nev %d -st_type sinvert -eps_gen_non_hermitian -st_pc_type lu -st_pc_factor_mat_solver_type superlu_dist", solver.nev)
    end

    #cool, this works, args are not picked up when we run julia -e ' ' arg1 etc though!
    for i in ARGS
        slepcargs *= " " * i
    end

    #hermitian may just not be true for qfm (or at all!). Perhaps this was causing some problemos?
    #this is at least effecting the inner product error a bit, unsure if it will fix it fully tbh
    #should always be hermitian, at least to tolerance, need to understand why not.
    #slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_non_hermitian -eps_view", solver.nev) #* evals_str #* efuncs_str 
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

    #fkn awful solution
    #think we need a separate inst_problem function we can call
    #removed, not going to work!
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
