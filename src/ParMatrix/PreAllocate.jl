
"""
Function that preallocates the number of non-zero elements in the petsc sparse matrices. This leads to memory efficient matrix construction. First the matrix is divided between processors by splitting the radial grid up. Then the number of non-zeros is computed per row, and this is split into the diagonal and the non-diagonal part, to conform with petsc's preallocation method.
"""
function preallocate_matrix(grids::MID.GridsT)


    local_n = split_matrix(grids)


    W = MatCreate()
    I = MatCreate()

    global_n = matrix_size(grids) 

    #this step is probably not needed, but just certifies that we know exactly which part each proc is
    #handling.
    MatSetSizes(W, local_n, local_n, global_n, global_n)
    MatSetSizes(I, local_n, local_n, global_n, global_n)

    MatSetFromOptions(W)
    MatSetFromOptions(I)

    indstart, indend = MatGetOwnershipRange(W)

    #@printf("Proc %d has %d to %d, local_n is %d\n", MPI.Comm_rank(MPI.COMM_WORLD), indstart, indend, local_n)



    inds = indstart:indend -1

    

    dnnz = zeros(Int32, local_n)
    onnz = zeros(Int32, local_n)


    #probbaly going to ignore the boundaries, will be much to complicated,
    #and for large grids should become negligible.
    boundary_inds = compute_boundary_inds(grids)

    
    

    #iterate through each row owned by this processor.
    for i in 1:local_n

        #plus one to change to julia form.
        if inds[i]+1 in boundary_inds
            #this only removes the boundaries from the rows not the columns
            #other boundaries are done inside compute_nz_inds, but this 
            #completly ignores rows, improving efficiency.
            dnnz[i] = 1
            onnz[i] = 0
            continue
        end 
        
        #determines the non-zero indicies for this row.
        nz_inds = compute_nz_inds(i, grids, inds, boundary_inds)
            
        

        #finds the number of nz's inside the processors diagonal,
        #see https://petsc.org/release/manualpages/Mat/MatMPIAIJSetPreallocation/
        #for example of diagonal and off-diagonal.
        dnnz[i] = length(nz_inds[indstart .< nz_inds .<= indend])
        onnz[i] = length(nz_inds[indstart .>= nz_inds])
        onnz[i] += length(nz_inds[nz_inds .> indend])





    end

    #allocates the memory for hte matrices.
    #MatMPIAIJSetPreallocation(W, PetscInt(1), dnnz, PetscInt(1), onnz)
    #MatMPIAIJSetPreallocation(I, PetscInt(1), dnnz, PetscInt(1), onnz)

    #finish the matrix setup.
    MatSetUp(W)
    MatSetUp(I)


    return W, I
end


#split the matrix up between the cores, setup such that each split is on a grid point.
#think we can only split on a fem grid point tbh.
function split_matrix(grids::MID.FFFGridsT)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm) #rank of each worker
    nprocs = MPI.Comm_size(comm) #total number of workers including root.
    root = 0
    counts = zeros(Int64, nprocs)
    #splits = zeros(Int64, nprocs+1)

    #here we split up the radial grid for each proc
    if rank==root
        grid_size = grids.r.N * grids.θ.N * grids.ζ.N
        counts_guess = Int64(div(grid_size, nprocs, RoundDown))
        Remainder    = Int64(grid_size - counts_guess*nprocs)
        counts[:]   .= counts_guess
        for i in 1:Remainder
            counts[i] += 1
        end
        #splits[1:end-1] = cumsum(append!([0], counts))[1:nprocs] 
        #splits[end] = grids.r.N - 1
    end

    #MPI.Bcast!(splits, root, comm)
    MPI.Bcast!(counts, root, comm)

    

    #block size is how many points indicies per radial point.
    local_n = counts[rank+1] * 8 # * grids.θ.N * grids.ζ.N

    return local_n
    #compute_block_size(grids)
end

function split_matrix(grids::MID.FFSGridsT)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm) #rank of each worker
    nprocs = MPI.Comm_size(comm) #total number of workers including root.
    root = 0
    counts = zeros(Int64, nprocs)
    #splits = zeros(Int64, nprocs+1)

    #here we split up the radial grid for each proc
    if rank==root
        grid_size = grids.r.N * grids.θ.N 
        counts_guess = Int64(div(grid_size, nprocs, RoundDown))
        Remainder    = Int64(grid_size - counts_guess*nprocs)
        counts[:]   .= counts_guess
        for i in 1:Remainder
            counts[i] += 1
        end
        #splits[1:end-1] = cumsum(append!([0], counts))[1:nprocs] 
        #splits[end] = grids.r.N - 1
    end

    #MPI.Bcast!(splits, root, comm)
    MPI.Bcast!(counts, root, comm)

    

    #block size is how many points indicies per radial point.
    local_n = counts[rank+1] * 4 * grids.ζ.N

    return local_n
    #compute_block_size(grids)
end

function split_matrix(grids::MID.FSSGridsT)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm) #rank of each worker
    nprocs = MPI.Comm_size(comm) #total number of workers including root.
    root = 0
    counts = zeros(Int64, nprocs)
    #splits = zeros(Int64, nprocs+1)

    #here we split up the radial grid for each proc
    if rank==root
        grid_size = grids.r.N 
        counts_guess = Int64(div(grid_size, nprocs, RoundDown))
        Remainder    = Int64(grid_size - counts_guess*nprocs)
        counts[:]   .= counts_guess
        for i in 1:Remainder
            counts[i] += 1
        end
        #splits[1:end-1] = cumsum(append!([0], counts))[1:nprocs] 
        #splits[end] = grids.r.N - 1
    end

    #MPI.Bcast!(splits, root, comm)
    MPI.Bcast!(counts, root, comm)

    

    #block size is how many points indicies per radial point.
    local_n = counts[rank+1] * 2 * grids.θ.N * grids.ζ.N

    return local_n
    #compute_block_size(grids)
end


