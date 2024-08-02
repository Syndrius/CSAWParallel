
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

include("PreAllocate.jl")



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
    slepcargs = @sprintf("-eps_nev %d -eps_target %s -st_type sinvert -memory_view -mat_view ::ascii_info", nev, target_freq) * evals_str * efuncs_str
    
    #this should also init petsc
    SlepcInitialize(slepcargs)

    #W, I = init_petsc_matrix(grids)
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display("Preparing Matrices...")
    end
    W, I = preallocate_matrix(grids)


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


#now works liks the other one, little be extra, but seems pretty good.
#this one is pretty good, as matrix is always pretty diagonal, 
#others are completly cooked lol.
function init_petsc_matrix(grids::MID.FSSGridsT)

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
    local_n = counts[rank+1] * grids.θ.count * grids.ζ.count * 2


    W = MatCreate()
    I = MatCreate()

    #this step is probably not needed, but just certifies that we know exactly which part each proc is
    #handling.
    MatSetSizes(W, local_n, local_n, global_n, global_n)
    MatSetSizes(I, local_n, local_n, global_n, global_n)

    MatSetFromOptions(W)
    MatSetFromOptions(I)

    indstart, indend = MatGetOwnershipRange(W)

    #indstart = indstart + 1 #get to julia form. indend does not need to change.
    

    display((indstart, indend))
    #plus one to shift to julia indexing
    #r_start = Int64(indstart / (grids.θ.N * grids.ζ.count * 4)) + 1
    #this has a plus 1 for indexing, and a minus 1 as ownership range returns end+1.
    #r_end = Int64(indend / (grids.θ.N * grids.ζ.N * 4))+1

    #rlocal = r_start:r_end

    indslocal = indstart:indend 


    #local diagonal matrix size.


    dnnz = zeros(Int32, local_n)
    onnz = zeros(Int32, local_n)

    #bc_edge1 = 8 * grids.θ.N * grids.ζ.N

    

    #bc_edge2 = global_n
    block_size = 2 * grids.θ.count * grids.ζ.count

    boundary_inds = MID.compute_boundary_inds(grids)


    for i in 1:local_n

        local_block = Int64(div(indslocal[i], block_size, RoundDown)) + 1

        #if rank==root
        #    display(local_block)
        #end

        #ie set the boundary cases to zero.
        #plus one needed to match julia indexing.
        if indslocal[i]+1 in boundary_inds
            dnnz[i] = 1
            onnz[i] = 0
            continue
        end

        if indslocal[i] < block_size 

            #this is like the opposite of the boundary inds
            boundary1 = 2:2:block_size
            #boundary2 = 4:4:block_size
            #boundary3 = 7:8:block_size
            #boundary4 = 8:8:block_size

            boundary2 = block_size+1:2*block_size

            #these are the non-zeros, avoiding the boundary.
            #stupid af name holy moley.
            nz_inds = vcat(boundary1, boundary2)


            #seems like one of these should be an >=??
            #weird we don't get an error.
            dnz = length(nz_inds[nz_inds .<= indend])
            onz = length(nz_inds[nz_inds .> indend])

        
        elseif indslocal[i] < 2 * block_size

            boundary1 = 2:2:block_size
            #boundary2 = 4:4:block_size
            #boundary3 = 7:8:block_size
            #boundary4 = 8:8:block_size

            boundary2 = block_size+1:3*block_size

            #these are the non-zeros, avoiding the boundary.
            #stupid af name holy moley.
            nz_inds = vcat(boundary1, boundary2)
            

            dnz = length(nz_inds[indstart .< nz_inds .<= indend])
            onz = length(nz_inds[indstart .>= nz_inds])
            onz += length(nz_inds[nz_inds .> indend])

           

        elseif indslocal[i] > global_n - block_size

            boundary1 = 2 + (grids.r.N - 1) * block_size:2:grids.r.N * block_size
            #boundary2 = 4 + (grids.r.N - 1) * block_size:4:grids.r.N * block_size
            #boundary3 = 7 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size
            #boundary4 = 8 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size

            boundary2 = (grids.r.N - 2) * block_size+1: (grids.r.N - 1) * block_size


            nz_inds = vcat(boundary1, boundary2)

            dnz = length(nz_inds[indstart .< nz_inds])
            onz = length(nz_inds[nz_inds .<= indstart])

            #nzs = Int32(block_size / 2) + block_size

            #nzs = 2 * block_size

            #nzstart = (local_block - 2) * block_size
            #nzend = (local_block+1) * block_size

            #nzstart = indslocal[i] - block_size
            #nzend = indslocal[i] + 2 * block_size

            #onz = max(nzend - indend, 0)
            #onz = max(indstart - nzstart, 0)
            #dnz = nzs - onz

        elseif indslocal[i] > global_n - 2* block_size

            boundary1 = 2 + (grids.r.N - 1) * block_size:2:grids.r.N * block_size
            #boundary2 = 4 + (grids.r.N - 1) * block_size:4:grids.r.N * block_size
            #boundary3 = 7 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size
            #boundary4 = 8 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size

            boundary2 = (grids.r.N - 3) * block_size+1: (grids.r.N - 1) * block_size


            nz_inds = vcat(boundary1, boundary2)

            dnz = length(nz_inds[indstart .< nz_inds .<= indend])
            onz = length(nz_inds[nz_inds .> indend])
            onz += length(nz_inds[nz_inds .<= indstart])


            #nzs = Int32(block_size / 2) + 2*block_size

            #nzs = 3 * block_size

            #nzstart = (local_block - 2) * block_size
            #nzend = (local_block+1) * block_size

            #nzstart = indslocal[i] - block_size
            #nzend = indslocal[i] + 2 * block_size

            #onz = max(nzend - indend, 0)
            #onz += max(indstart - nzstart, 0)
            #dnz = nzs - onz

        else
            nzs = 3 * block_size

            nzstart = (local_block - 2) * block_size
            nzend = (local_block+1) * block_size

            #nzstart = indslocal[i] - block_size
            #nzend = indslocal[i] + 2 * block_size

            onz = max(nzend - indend, 0)
            onz += max(indstart - nzstart, 0)
            dnz = nzs - onz

        end

        #this essentially ignores the boundary conditions.
        

        dnnz[i] = dnz
        onnz[i] = onz
    end

    if rank==1
        println(dnnz)
        println(onnz)
    end

    
    MatMPIAIJSetPreallocation(W, PetscInt(1), dnnz, PetscInt(1), onnz)
    MatMPIAIJSetPreallocation(I, PetscInt(1), dnnz, PetscInt(1), onnz)


    MatSetUp(W)
    MatSetUp(I)


    return W, I

end


#now works liks the other one, little be extra, but seems pretty good.
#seems to require θ.count * ζ.count extra for some reason, but pretty much perfect.
function init_petsc_matrix(grids::MID.FFSGridsT)

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
    local_n = counts[rank+1] * grids.θ.N * grids.ζ.count * 4


    W = MatCreate()
    I = MatCreate()

    #this step is probably not needed, but just certifies that we know exactly which part each proc is
    #handling.
    MatSetSizes(W, local_n, local_n, global_n, global_n)
    MatSetSizes(I, local_n, local_n, global_n, global_n)

    MatSetFromOptions(W)
    MatSetFromOptions(I)

    indstart, indend = MatGetOwnershipRange(W)

    #indstart = indstart + 1 #get to julia form. indend does not need to change.
    

    
    #plus one to shift to julia indexing
    #r_start = Int64(indstart / (grids.θ.N * grids.ζ.count * 4)) + 1
    #this has a plus 1 for indexing, and a minus 1 as ownership range returns end+1.
    #r_end = Int64(indend / (grids.θ.N * grids.ζ.N * 4))+1

    #rlocal = r_start:r_end

    indslocal = indstart:indend 

    #local diagonal matrix size.


    dnnz = zeros(Int32, local_n)
    onnz = zeros(Int32, local_n)

    #bc_edge1 = 8 * grids.θ.N * grids.ζ.N

    

    #bc_edge2 = global_n
    block_size = 4 * grids.θ.N * grids.ζ.count

    boundary_inds = MID.compute_boundary_inds(grids)


    for i in 1:local_n

        local_block = Int64(div(indslocal[i], block_size, RoundDown)) + 1

        #if rank==root
        #    display(local_block)
        #end

        #ie set the boundary cases to zero.
        #plus one needed to match julia indexing.
        if indslocal[i]+1 in boundary_inds
            dnnz[i] = 1
            onnz[i] = 0
            continue
        end

        if indslocal[i] < block_size 

            #this is like the opposite of the boundary inds
            boundary1 = 3:4:block_size
            boundary2 = 4:4:block_size
            #boundary3 = 7:8:block_size
            #boundary4 = 8:8:block_size

            boundary3 = block_size+1:2*block_size

            #these are the non-zeros, avoiding the boundary.
            #stupid af name holy moley.
            nz_inds = vcat(boundary1, boundary2, boundary3)


            #seems like one of these should be an >=??
            #weird we don't get an error.
            dnz = length(nz_inds[nz_inds .<= indend])
            onz = length(nz_inds[nz_inds .> indend])

        
        elseif indslocal[i] < 2 * block_size

            boundary1 = 3:4:block_size
            boundary2 = 4:4:block_size
            #boundary3 = 7:8:block_size
            #boundary4 = 8:8:block_size

            boundary3 = block_size+1:3*block_size

            #these are the non-zeros, avoiding the boundary.
            #stupid af name holy moley.
            nz_inds = vcat(boundary1, boundary2, boundary3)
            

            dnz = length(nz_inds[indstart .< nz_inds .<= indend])
            onz = length(nz_inds[indstart .>= nz_inds])
            onz += length(nz_inds[indend .< nz_inds])

           

        elseif indslocal[i] > global_n - block_size

            boundary1 = 3 + (grids.r.N - 1) * block_size:4:grids.r.N * block_size
            boundary2 = 4 + (grids.r.N - 1) * block_size:4:grids.r.N * block_size
            #boundary3 = 7 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size
            #boundary4 = 8 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size

            boundary3 = (grids.r.N - 2) * block_size+1: (grids.r.N - 1) * block_size


            nz_inds = vcat(boundary1, boundary2, boundary3)

            dnz = length(nz_inds[indstart .< nz_inds])
            onz = length(nz_inds[nz_inds .<= indstart])

            #nzs = Int32(block_size / 2) + block_size

            #nzs = 2 * block_size

            #nzstart = (local_block - 2) * block_size
            #nzend = (local_block+1) * block_size

            #nzstart = indslocal[i] - block_size
            #nzend = indslocal[i] + 2 * block_size

            #onz = max(nzend - indend, 0)
            #onz = max(indstart - nzstart, 0)
            #dnz = nzs - onz

        elseif indslocal[i] > global_n - 2* block_size

            boundary1 = 3 + (grids.r.N - 1) * block_size:4:grids.r.N * block_size
            boundary2 = 4 + (grids.r.N - 1) * block_size:4:grids.r.N * block_size
            #boundary3 = 7 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size
            #boundary4 = 8 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size

            boundary3 = (grids.r.N - 3) * block_size+1: (grids.r.N - 1) * block_size


            nz_inds = vcat(boundary1, boundary2, boundary3)

            dnz = length(nz_inds[indstart .< nz_inds .<= indend])
            onz = length(nz_inds[nz_inds .> indend])
            onz += length(nz_inds[nz_inds .<= indstart])


            #nzs = Int32(block_size / 2) + 2*block_size

            #nzs = 3 * block_size

            #nzstart = (local_block - 2) * block_size
            #nzend = (local_block+1) * block_size

            #nzstart = indslocal[i] - block_size
            #nzend = indslocal[i] + 2 * block_size

            #onz = max(nzend - indend, 0)
            #onz += max(indstart - nzstart, 0)
            #dnz = nzs - onz

        else
            nzs = 3 * block_size

            nzstart = (local_block - 2) * block_size
            nzend = (local_block+1) * block_size

            #nzstart = indslocal[i] - block_size
            #nzend = indslocal[i] + 2 * block_size

            onz = max(nzend - indend, 0)
            onz += max(indstart - nzstart, 0)
            dnz = nzs - onz

        end

        #this essentially ignores the boundary conditions.
        

        dnnz[i] = dnz
        onnz[i] = onz
    end


    
    MatMPIAIJSetPreallocation(W, PetscInt(1), dnnz, PetscInt(1), onnz)
    MatMPIAIJSetPreallocation(I, PetscInt(1), dnnz, PetscInt(1), onnz)


    MatSetUp(W)
    MatSetUp(I)


    return W, I

end


#only works for FFF, ideally will have one of these for each method.
function init_petsc_matrix(grids::MID.FFFGridsT)

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

    #Need to compute the non-zero's per row properly.

    #these define the diagonal and off-diagonal non-zeros for each row.
    



    

    W = MatCreate()
    I = MatCreate()

    #this step is probably not needed, but just certifies that we know exactly which part each proc is
    #handling.
    MatSetSizes(W, local_n, local_n, global_n, global_n)
    MatSetSizes(I, local_n, local_n, global_n, global_n)

    MatSetFromOptions(W)
    MatSetFromOptions(I)

    indstart, indend = MatGetOwnershipRange(W)
    #plus one to shift to julia indexing
    #r_start = Int64(indstart / (grids.θ.N * grids.ζ.N * 8)) + 1
    #this has a plus 1 for indexing, and a minus 1 as ownership range returns end+1.
    #r_end = Int64(indend / (grids.θ.N * grids.ζ.N * 8))+1

    #rlocal = r_start:r_end

    indslocal = indstart:indend 

    #local diagonal matrix size.


    dnnz = zeros(Int32, local_n)
    onnz = zeros(Int32, local_n)

    #bc_edge1 = 8 * grids.θ.N * grids.ζ.N

    

    #bc_edge2 = global_n
    block_size = 8 * grids.θ.N * grids.ζ.N

    boundary_inds = MID.compute_boundary_inds(grids)

    #non_boundary_inds = deleteat!(collect(1:block_size), sort(unique(boundary_inds)))

    
    #if rank==0
        #display(boundary_inds)
        #display(nond_boundary_inds)
    #end
    #this is working well enough now, still predicting a few extra values, 
    #but I think this is just extra zeros.
    for i in 1:local_n

        local_block = Int64(div(indslocal[i], block_size, RoundDown)) + 1

        #display(local_block)

        #ie set the boundary cases to zero.
        if indslocal[i]+1 in boundary_inds
            dnnz[i] = 1
            onnz[i] = 0
            continue
        end

        if indslocal[i] < block_size 

            #this is like the opposite of the boundary inds
            boundary1 = 5:8:block_size
            boundary2 = 6:8:block_size
            boundary3 = 7:8:block_size
            boundary4 = 8:8:block_size

            boundary5 = block_size+1:2*block_size

            #these are the non-zeros, avoiding the boundary.
            #stupid af name holy moley.
            nz_inds = vcat(boundary1, boundary2, boundary3, boundary4, boundary5)

            #start by working out how many nz there will be, then we can decide if they are in the diagonal or not.
            #nzs = Int32(block_size / 2) + block_size

            #dnz = count(j->(j<indend), [non_boundary_inds])
            #onz = count(j->(j>indend), [non_boundary_inds])

            #seems like one of these should be an >=??
            #weird we don't get an error.
            dnz = length(nz_inds[nz_inds .<= indend])
            onz = length(nz_inds[nz_inds .> indend])

            #nzs = 2 * block_size

            #nzstart = 0#(local_block - 2) * block_size
            #nzend = (local_block+1) * block_size

            #nzstart = indslocal[i] - block_size
            #nzend = indslocal[i] + 2 * block_size

            #onz = max(nzend - indend, 0)
            #onz += max(indstart - nzstart, 0)
            #dnz = nzs - onz



            #do stuff
        
        elseif indslocal[i] < 2 * block_size

            boundary1 = 5:8:block_size
            boundary2 = 6:8:block_size
            boundary3 = 7:8:block_size
            boundary4 = 8:8:block_size

            boundary5 = block_size+1:3*block_size

            #these are the non-zeros, avoiding the boundary.
            #stupid af name holy moley.
            nz_inds = vcat(boundary1, boundary2, boundary3, boundary4, boundary5)

            dnz = length(nz_inds[indstart .< nz_inds .<= indend])
            onz = length(nz_inds[indstart .>= nz_inds])
            onz += length(nz_inds[nz_inds .> indend])

            #nzs = Int32(block_size / 2) + 2*block_size

            #half the normal amount of non-zeros
            #nzs = Int32(3 * block_size / 2)

            #nzs = 3 * block_size

            #nzstart = (local_block - 2) * block_size
            #nzend = (local_block+1) * block_size

            #nzstart = indslocal[i] - block_size
            #nzend = indslocal[i] + 2 * block_size

            #onz = max(nzend - indend, 0)
            #onz += max(indstart - nzstart, 0)
            #dnz = nzs - onz

        elseif indslocal[i] > global_n - block_size

            boundary1 = 5 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size
            boundary2 = 6 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size
            boundary3 = 7 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size
            boundary4 = 8 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size

            boundary5 = (grids.r.N - 2) * block_size+1: (grids.r.N - 1) * block_size


            nz_inds = vcat(boundary1, boundary2, boundary3, boundary4, boundary5)

            dnz = length(nz_inds[nz_inds .> indstart])
            onz = length(nz_inds[nz_inds .<= indstart])

            #nzs = Int32(block_size / 2) + block_size

            #nzs = 2 * block_size

            #nzstart = (local_block - 2) * block_size
            #nzend = (local_block+1) * block_size

            #nzstart = indslocal[i] - block_size
            #nzend = indslocal[i] + 2 * block_size

            #onz = max(nzend - indend, 0)
            #onz = max(indstart - nzstart, 0)
            #dnz = nzs - onz

        elseif indslocal[i] > global_n - 2* block_size

            boundary1 = 5 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size
            boundary2 = 6 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size
            boundary3 = 7 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size
            boundary4 = 8 + (grids.r.N - 1) * block_size:8:grids.r.N * block_size

            boundary5 = (grids.r.N - 3) * block_size+1: (grids.r.N - 1) * block_size


            nz_inds = vcat(boundary1, boundary2, boundary3, boundary4, boundary5)

            dnz = length(nz_inds[indstart .< nz_inds .<= indend])
            onz = length(nz_inds[nz_inds .> indend])
            onz += length(nz_inds[nz_inds .<= indstart])


            #nzs = Int32(block_size / 2) + 2*block_size

            #nzs = 3 * block_size

            #nzstart = (local_block - 2) * block_size
            #nzend = (local_block+1) * block_size

            #nzstart = indslocal[i] - block_size
            #nzend = indslocal[i] + 2 * block_size

            #onz = max(nzend - indend, 0)
            #onz += max(indstart - nzstart, 0)
            #dnz = nzs - onz

        else
            nzs = 3 * block_size

            nzstart = (local_block - 2) * block_size
            nzend = (local_block+1) * block_size

            #nzstart = indslocal[i] - block_size
            #nzend = indslocal[i] + 2 * block_size

            onz = max(nzend - indend, 0)
            onz += max(indstart - nzstart, 0)
            dnz = nzs - onz

        end

        #this essentially ignores the boundary conditions.
        

        dnnz[i] = dnz
        onnz[i] = onz
    end

    #display(dnnz)
    #display(onnz)
            #need to check if we are at the edge of a local matrix now...

    #block is how many indicies per radial point.
    
    
    #each proc will have two non-zero blocks, one at the start, one at the end, we
    #can probably do this better, but this is a good start.

    #overestimate, only at boundaries between procs will there be any off-diagonal
    #and there it will be a single block (block_size per row for block_size columns).
    #nzo = block_size 

    #overestimate, at at boundaries between procs and true boundaries, there will typically only be 2 * block_size
    #we are also ignoring all the zero's introduced by boundary conditions.
    #nzd = 3 * block_size
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
    #significantly better
    #still need to deal with the boundary conditions though!
    MatMPIAIJSetPreallocation(W, PetscInt(1), dnnz, PetscInt(1), onnz)
    MatMPIAIJSetPreallocation(I, PetscInt(1), dnnz, PetscInt(1), onnz)
    #MatMPIAIJSetPreallocation(I, PetscInt(nzd), PetscInt(nzo))

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