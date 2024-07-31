
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

    #new version that deals with the petsc matrix better.

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
    slepcargs = @sprintf("-eps_nev %d -eps_target %s -st_type sinvert", nev, target_freq) * evals_str * efuncs_str
    
    #this should also init petsc
    SlepcInitialize(slepcargs)

    W, I = init_petsc_matrix(grids)


    #construct the matrices, each proc has a subset of the coo data.        
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Constructing...")
    end
    par_construct(W, I, prob, grids)

    MatAssemblyBegin(W, MAT_FINAL_ASSEMBLY)
    MatAssemblyBegin(I, MAT_FINAL_ASSEMBLY)
    MatAssemblyEnd(W, MAT_FINAL_ASSEMBLY)
    MatAssemblyEnd(I, MAT_FINAL_ASSEMBLY)

    

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Solving...")
    end
    #memory doubling does not seem to be in solve...
    #solve the matrix, this writes to file
    par_solve(W, I)

    destroy!(W)
    destroy!(I)
    SlepcFinalize()
    MPI.Finalize()

end


#only works for FFF, ideally will have one of these for each method.
function init_petsc_matrix(grids)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm) #rank of each worker
    nprocs = MPI.Comm_size(comm) #total number of workers including root.
    root = 0
    counts = zeros(Int64, nprocs)
    #splits = zeros(Int64, nprocs+1)

    #here we split up the radial grid for each proc
    if rank==root
    
        counts_guess = Int64(div(grids.r.N, nprocs, RoundDown))
        Remainder    = Int64(grids.r.N - counts_guess*nprocs)
        counts[:]   .= counts_guess
        for i in 1:Remainder
            counts[i] += 1
        end
        #splits[1:end-1] = cumsum(append!([0], counts))[1:nprocs] 
        #splits[end] = grids.r.N - 1
    end

    #MPI.Bcast!(splits, root, comm)
    MPI.Bcast!(counts, root, comm)

    #if rank==0
    #    display("counts")
    #    display(counts)
    #end
    #compute the matrix size for creating the Petsc matrix.
    global_n = matrix_dim(grids) #should work for both grids now!
    #n = 2 * grids.r.N * grids.pmd.count * grids.tmd.count
    #W = create_matrix(global_n, global_n, auto_setup=true)
    #I = create_matrix(global_n, global_n, auto_setup=true)
    #we will need a function for this, and it will depend on the type of grid
    #probably just stick with FFF for now.
    #rank starts from 0...
    local_n = counts[rank+1] * grids.θ.N * grids.ζ.N * 8
    #if rank==0
        #display("local")
        #display(local_n)
        #display("global")
        #display(global_n)
        
    #end
    #this returns start of local range, and end+1 of local range, not sure why.
    #display(MatGetOwnershipRange(W))

    

    W = MatCreate()
    I = MatCreate()

    #this step is probably not needed, but just certifies that we know exactly which part each proc is
    #handling.
    MatSetSizes(W, local_n, local_n, global_n, global_n)
    MatSetSizes(I, local_n, local_n, global_n, global_n)

    MatSetFromOptions(W)
    MatSetFromOptions(I)

    #block is how many indicies per radial point.
    block_size = 8 * grids.θ.N * grids.ζ.N
    
    #each proc will have two non-zero blocks, one at the start, one at the end, we
    #can probably do this better, but this is a good start.

    #overestimate, only at boundaries between procs will there be any off-diagonal
    #and there it will be a single block (block_size per row for block_size columns).
    nzo = block_size 

    #overestimate, at at boundaries between procs and true boundaries, there will typically only be 2 * block_size
    #we are also ignoring all the zero's introduced by boundary conditions.
    nzd = 3 * block_size
    #blocks_per_proc = Int64(local_n / block_size)
    

    #display("local n")
    #display(local_n)
    #if rank==0
    #    display("block size")
    #    display(block_size)
    #end

    #display(mod(local_n, block_size))

    #this works, adequatly, but causes problems when not enough processors are used, 
    #ie the local range overlaps the block size, May need to assert that the local range is 
    #smaller than the block size...
    #display((nzd, nzo))
    MatMPIAIJSetPreallocation(W, PetscInt(nzd), PetscInt(nzo))
    MatMPIAIJSetPreallocation(I, PetscInt(nzd), PetscInt(nzo))

    MatSetUp(W)
    MatSetUp(I)


    return W, I
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
#new version that uses petsc matrices entirely to try and reduce memory usage.
function no_pre_par_construct_and_solve(; prob::MID.ProblemT, grids::MID.GridsT, σ=0.0::Float64, nev=50::Int64, dir::String)

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
    slepcargs = @sprintf("-eps_nev %d -eps_target %s -st_type sinvert -eps_view", nev, target_freq) * evals_str * efuncs_str
    
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
    par_construct(W, I, prob=prob, grids=grids)


    

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Solving...")
    end
    #memory doubling does not seem to be in solve...
    #solve the matrix, this writes to file
    par_solve(W, I)

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