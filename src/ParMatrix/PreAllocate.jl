
"""
    preallocate_matrix(grids::GridsT)

Function that preallocates the number of non-zero elements in the petsc sparse matrices. This leads to memory efficient matrix construction. First the matrix is divided between processors by splitting the radial grid up. Then the number of non-zeros is computed per row, and this is split into the diagonal and the non-diagonal part, to conform with petsc's preallocation method.
"""
function preallocate_matrix(grids::GridsT)

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

    #needs to be moved for flr effects.
    PetscWrap.MatSetOption(W, PetscWrap.MAT_HERMITIAN, true)
    PetscWrap.MatSetOption(I, PetscWrap.MAT_HERMITIAN, true)

    indstart, indend = MatGetOwnershipRange(W)

    #makes the final ind inclusive
    inds = indstart:indend -1

    dnnz = zeros(Int32, local_n)
    onnz = zeros(Int32, local_n)

    #gets the boundary inds, as these will be zero and therefore not allocated.
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

    #allocates the memory for the matrices.
    MatMPIAIJSetPreallocation(W, PetscInt(1), dnnz, PetscInt(1), onnz)
    MatMPIAIJSetPreallocation(I, PetscInt(1), dnnz, PetscInt(1), onnz)

    #finish the matrix setup.
    MatSetUp(W)
    MatSetUp(I)

    return W, I
end


"""
    split_matrix(grids::FFFGridsT)

Divides the matrix between cores.
"""
function split_matrix(grids::FFFGridsT)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm) #rank of each worker
    nprocs = MPI.Comm_size(comm) #total number of workers including root.
    root = 0
    counts = zeros(Int64, nprocs)

    #here we split up the radial grid for each proc
    if rank==root
        grid_size = grids.x1.N * grids.x2.N * grids.x3.N
        counts_guess = Int64(div(grid_size, nprocs, RoundDown))
        Remainder    = Int64(grid_size - counts_guess*nprocs)
        counts[:]   .= counts_guess
        for i in 1:Remainder
            counts[i] += 1
        end
    end

    MPI.Bcast!(counts, root, comm)

    #block size is how many indicies per grid point.
    local_n = counts[rank+1] * 8 

    return local_n
end


"""
    split_matrix(grids::FFSGridsT)

Divides the matrix between cores.
"""
function split_matrix(grids::FFSGridsT)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm) #rank of each worker
    nprocs = MPI.Comm_size(comm) #total number of workers including root.
    root = 0
    counts = zeros(Int64, nprocs)

    #here we split up the radial grid for each proc
    if rank==root
        grid_size = grids.x1.N * grids.x2.N 
        counts_guess = Int64(div(grid_size, nprocs, RoundDown))
        Remainder    = Int64(grid_size - counts_guess*nprocs)
        counts[:]   .= counts_guess
        for i in 1:Remainder
            counts[i] += 1
        end
    end

    MPI.Bcast!(counts, root, comm)

    #block size is how many indicies per grid point.
    local_n = counts[rank+1] * 4 * grids.x3.N

    return local_n
end


"""
    split_matrix(grids::FSSGridsT)

Divides the matrix between cores.
"""
function split_matrix(grids::FSSGridsT)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm) #rank of each worker
    nprocs = MPI.Comm_size(comm) #total number of workers including root.
    root = 0
    counts = zeros(Int64, nprocs)

    #here we split up the radial grid for each proc
    if rank==root
        grid_size = grids.x1.N 
        counts_guess = Int64(div(grid_size, nprocs, RoundDown))
        Remainder    = Int64(grid_size - counts_guess*nprocs)
        counts[:]   .= counts_guess
        for i in 1:Remainder
            counts[i] += 1
        end
    end

    MPI.Bcast!(counts, root, comm)

    #block size is how many indicies per grid point.
    local_n = counts[rank+1] * 2 * grids.x2.N * grids.x3.N

    return local_n
end
