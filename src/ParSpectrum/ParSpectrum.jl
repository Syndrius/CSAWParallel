
module ParSpectrum

#so this works, obvs without registering MID, we have to add this locally, ie Pkg.add ~/phd/MID or whatever
using MID #must be added locally as we do not have MID on the registary yet.
using MID.Indexing
using MID.Basis
using MID.WeakForm
using MID.Integration
using MID.Structures
using MPI
using FFTW
using FastGaussQuadrature
using SparseArrays
using Printf
using JLD2 #allows storage of 3d matrix in simple julia format.


#note to get these files to work on mac, had to modify the load.jl files in both cases
#mac uses dylib files, while linux uses .so files so PetscWrap and SlepcWrap were unable to find the files.
using PetscWrap
using SlepcWrap


export par_compute_spectrum
export par_spectrum_from_file


include("ParMatrixToGrid.jl")

include("ParConstruct.jl")


include("ParSolve.jl")


include("PreAllocate.jl")


include("ParPostProcess.jl")

export process_hdf5 #ideally this would not be needed but some weird shit is happening.
export process_hdf5_deriv



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
function par_compute_spectrum(; prob::MID.ProblemT, grids::MID.GridsT, target_freq=0.0::Float64, nev=200::Int64, dir::String)

    #un-normalise the target frequency for the shift and invert
    target_freq = target_freq^2 / prob.geo.R0^2 


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
    slepcargs = @sprintf("-eps_nev %d -eps_target %s -st_type sinvert -memory_view -mat_view ::ascii_info", nev, target_freq) * evals_str #* efuncs_str 

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
        mat_dim = MID.matrix_dim(grids)
        @printf("Construction of %dx%d matrices complete.\n", mat_dim, mat_dim)
    end


    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Solving...")
    end

    #solve the generalised eigenvalue problem with W and I.
    eps = par_solve(W, I)

    nconv = EPSGetConverged(eps)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @printf("Solving complete, %d eigenvalues found.\n", nconv)
    end

    #so this doesn't always work, may still be worth using, as it does seem to work most of the tim...
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Postprocessing...")
    end
    
    #vectors are used to store the eigenfunctions
    vecr, veci = MatCreateVecs(W)

    par_post_process(eps, dir, vecr, veci, nconv)

    #process the evals and efuncs into practical formats.
    #old_par_post_process(eps, dir, grids, prob.geo, vecr, veci, nconv)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Finished.")
    end

    #free memory and finalise.
    destroy!(vecr)
    destroy!(veci)
    destroy!(eps)
    destroy!(W)
    destroy!(I)
    SlepcFinalize()
    MPI.Finalize()

end


"""
    par_spectrum_from_file(; dir::String, target_freq::Float64, nev=100::Int64)

Computes the spectrum in parallel from inputs stored in files in the given directory.
"""
function par_spectrum_from_file(; dir::String, target_freq::Float64, nev=100::Int64)

    #should only root be doing this?
    prob, grids = MID.inputs_from_file(dir=dir)

    par_compute_spectrum(prob=prob, grids=grids, target_freq=target_freq, nev=nev, dir=dir)

end

end