
module ParSolve

#so this works, obvs without registering MID, we have to add this locally, ie Pkg.add ~/phd/MID or whatever
using Printf
using JLD2 #allows storage of 3d matrix in simple julia format.


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
export par_spectrum_from_file
export qfm_spectrum_from_file




include("PostProcess.jl")

export par_post_process #ideally this would not be needed but some weird shit is happening.
#export process_hdf5_deriv

include("ShiftInvertSolve.jl")


include("SliceSolve.jl")

export par_solve


"""
    par_compute_spectrum(; prob::MID.ProblemT, grids::MID.GridsT, target_freq=0.0::Float64, nev=100::Int64, dir::String)

Computes the spectrum in parallel, and writes solution to file.

# Args
- prob::ProblemT Struct containing the functions and parameters that define the problem we are solving
- grids::GridT Grids to solve over.
- target_freq::Float64 Frequency for 'shift and invert', typically this is the TAE frequency.
- nev=200::Int64 Number of eigenvalues to solve for.
- dir::String Directory the results are written to.
"""
function par_compute_spectrum(; prob::ProblemT, grids::GridsT, solver::SolverT, dir::String)

    #un-normalise the target frequency for the shift and invert
    #target_freq = target_freq^2 / prob.geo.R0^2 


    MPI.Init()
    #won't be writing the efuncs like this, as it is cooked.
    #efuncs_str = " -eps_view_vectors :" * dir * "funcs.dat:ascii_symmodu"
    evals_str = " -eps_view_values :" * dir * "vals.dat:ascii_matlab"
    #efuncs_str = " -eps_view_vectors ascii_python:" * dir * "funcs"
    #alternative args
    #eps_harmonics means we are searching inside the spectrum.
    #-st_type sinvert
    #slepcargs = @sprintf("-eps_nev %d -eps_target %s -eps_harmonic", nev, σ) * evals_str * efuncs_str
    #-eps_view for solver stuff
    #-memory_view for mem
    #log_view for heaps of petsc info.
    #need to make it auto detect if hermitian or not
    #slepcargs = @sprintf("-eps_nev %d -eps_target %s -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_non_hermitian -eps_view", nev, target_freq) * evals_str #* efuncs_str 

    #probably want this to consider the solver a bit, but probably doesn't matter for now.
    #don't think the solver is actually different for Hermitian cases.
    #removed eps target from here, so that that is set later!
    slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_non_hermitian -eps_view", solver.nev) * evals_str #* efuncs_str 



    ############
    #attempt at using CISS (https://slepc.upv.es/documentation/reports/str11.pdf)
    #for finding all evals in range rather than a specific target.
    #freq_low = 0.3^2 / prob.geo.R0^2
    #freq_low = -10
    #freq_high = 12^2 / prob.geo.R0^2
    #slepcargs = @sprintf("-eps_type ciss -rg_type interval -rg_interval_endpoints %f,%f,-10,10 -eps_ciss_maxblocksize 4000 -eps_ciss_integration_points 1000", freq_low, freq_high) * evals_str

    #slepcargs = @sprintf("-eps_interval %s,%s -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky -st_pc_factor_mat_solver_type superlu_dist -st_mat_superlu_dist_rowperm NOROWPERM", freq_low, freq_high) * evals_str
    
    #slepcargs = @sprintf("-eps_type ciss -st_type sinvert -rg_type interval -rg_interval_endpoints -10,10,-10,10") * evals_str
    #display(slepcargs)
    #slepcargs = @sprintf("-eps_type ciss -rg_type ellipse -rg_ellipse_center %s -mat_view ::ascii_info", target_freq) * evals_str


    #so this will work, we just have to be realistic about the integration range,
    #i guess it does use discrete points so would skip over many values if we consider an enourmous region.
    #this will take some work to be of practical use.
    #however may be w better option in the long run. Hard to know these things.
    #slepcargs = @sprintf("-eps_type ciss -rg_type ellipse -rg_ellipse_center 0.0014 -rg_ellipse_radius 0.002 -rg_ellipse_vscale 1.0 -st_pc_factor_shift_type NONZERO") * evals_str



    

    #may need to try this outside of julia... with debug on, no fkn idea what is going on.

    #initialise slepc, setting the number of eigenvalues (nev) the target frequency and declaring that shift and inver should be used.
    #slepcargs = @sprintf("-eps_nev %d -eps_target %s -st_type sinvert -mat_view ::ascii_info", nev, target_freq)
    #this should also init petsc
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
    #destroy!(eps)
    destroy!(W)
    destroy!(I)
    SlepcFinalize()
    MPI.Finalize()

end


#this function needs to not be an independent function!
#this function is probably telling us that the surfs belong in the problem!
function qfm_compute_spectrum(; prob::ProblemT, grids::GridsT, solver::SolverT, dir::String, surfs::Array{QFMSurfaceT})



    MPI.Init()
    #won't be writing the efuncs like this, as it is cooked.
    #efuncs_str = " -eps_view_vectors :" * dir * "funcs.dat:ascii_symmodu"
    #evals_str = " -eps_view_values :" * dir * "vals.dat:ascii_matlab"
    #efuncs_str = " -eps_view_vectors ascii_python:" * dir * "funcs"
    #alternative args
    #eps_harmonics means we are searching inside the spectrum.
    #-st_type sinvert
    #slepcargs = @sprintf("-eps_nev %d -eps_target %s -eps_harmonic", nev, σ) * evals_str * efuncs_str
    #-eps_view for solver stuff
    #-memory_view for mem
    #log_view for heaps of petsc info.
    #need to make it auto detect if hermitian or not
    #probably want to actually set nev later. Maybe best if much of this is not actually set here manually. rather it should be set based on the solver type
    #eg non-hermitain etc.
    slepcargs = @sprintf("-eps_nev %d -st_type sinvert -memory_view -mat_view ::ascii_info -eps_gen_non_hermitian -eps_view", solver.nev) #* evals_str #* efuncs_str 


    ############
    #attempt at using CISS (https://slepc.upv.es/documentation/reports/str11.pdf)
    #for finding all evals in range rather than a specific target.
    #freq_low = 0.3^2 / prob.geo.R0^2
    #freq_low = -10
    #freq_high = 12^2 / prob.geo.R0^2
    #slepcargs = @sprintf("-eps_type ciss -rg_type interval -rg_interval_endpoints %f,%f,-10,10 -eps_ciss_maxblocksize 4000 -eps_ciss_integration_points 1000", freq_low, freq_high) * evals_str

    #slepcargs = @sprintf("-eps_interval %s,%s -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky -st_pc_factor_mat_solver_type superlu_dist -st_mat_superlu_dist_rowperm NOROWPERM", freq_low, freq_high) * evals_str
    
    #slepcargs = @sprintf("-eps_type ciss -st_type sinvert -rg_type interval -rg_interval_endpoints -10,10,-10,10") * evals_str
    #display(slepcargs)
    #slepcargs = @sprintf("-eps_type ciss -rg_type ellipse -rg_ellipse_center %s -mat_view ::ascii_info", target_freq) * evals_str


    #so this will work, we just have to be realistic about the integration range,
    #i guess it does use discrete points so would skip over many values if we consider an enourmous region.
    #this will take some work to be of practical use.
    #however may be w better option in the long run. Hard to know these things.
    #slepcargs = @sprintf("-eps_type ciss -rg_type ellipse -rg_ellipse_center 0.0014 -rg_ellipse_radius 0.002 -rg_ellipse_vscale 1.0 -st_pc_factor_shift_type NONZERO") * evals_str



    

    #may need to try this outside of julia... with debug on, no fkn idea what is going on.

    #initialise slepc, setting the number of eigenvalues (nev) the target frequency and declaring that shift and inver should be used.
    #slepcargs = @sprintf("-eps_nev %d -eps_target %s -st_type sinvert -mat_view ::ascii_info", nev, target_freq)
    #this should also init petsc
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
    #destroy!(eps)
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
