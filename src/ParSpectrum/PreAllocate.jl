
"""
Function that preallocates the number of non-zero elements in the petsc sparse matrices. This leads to memory efficient matrix construction. First the matrix is divided between processors by splitting the radial grid up. Then the number of non-zeros is computed per row, and this is split into the diagonal and the non-diagonal part, to conform with petsc's preallocation method.
"""
function preallocate_matrix(grids::MID.GridsT)


    local_n = split_matrix(grids)


    W = MatCreate()
    I = MatCreate()

    global_n = matrix_dim(grids) 

    #this step is probably not needed, but just certifies that we know exactly which part each proc is
    #handling.
    MatSetSizes(W, local_n, local_n, global_n, global_n)
    MatSetSizes(I, local_n, local_n, global_n, global_n)

    MatSetFromOptions(W)
    MatSetFromOptions(I)

    indstart, indend = MatGetOwnershipRange(W)

    #@printf("Proc %d has %d to %d, local_n is %d\n", MPI.Comm_rank(MPI.COMM_WORLD), indstart, indend, local_n)



    inds = indstart:indend 



    dnnz = zeros(Int32, local_n)
    onnz = zeros(Int32, local_n)


    #probbaly going to ignore the boundaries, will be much to complicated,
    #and for large grids should become negligible.
    boundary_inds = MID.compute_boundary_inds(grids)

    #if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        #display(boundary_inds)
    #end

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
    MatMPIAIJSetPreallocation(W, PetscInt(1), dnnz, PetscInt(1), onnz)
    MatMPIAIJSetPreallocation(I, PetscInt(1), dnnz, PetscInt(1), onnz)

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
    local_n = counts[rank+1] * 4 * grids.ζ.count

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
    local_n = counts[rank+1] * 2 * grids.θ.count * grids.ζ.count

    return local_n
    #compute_block_size(grids)
end


"""

Computes the indicies of the non-zero elements for a given row of the matrix. The non-zeros are split int `blocks` for each radial point. The blocks follow a diagonal pattern. In FSS case, each block is dense.
"""
function compute_nz_inds(ind, grids::MID.FSSGridsT, inds, boundary_inds)


    #size of the largest block in the sparse matrix,
    #essentially the numebr of points for each r grid point.
    block_size = compute_block_size(grids)

    #which row of blocks we are currently in.
    block_row = div(inds[ind], block_size) + 1

    #total number of blocks, which is the numebr of radial grids.
    max_block_row = grids.r.N

    #for fss case, each block is dense and fully occupied.
    #this is kept here to maintain consistency with other cases.
    sub_nz_inds = 1:block_size

    #If we are in the first row, we must deal with the boundary conditions,
    #and there will only be 2 blocks.
    if block_row == 1

        nz_inds1 = filter(x->  !(x in boundary_inds),  sub_nz_inds)
        nz_inds2 = block_size .+ sub_nz_inds
        nz_inds = vcat(nz_inds1, nz_inds2)

    #second row also deals with boundary conditions
    elseif block_row == 2
        nz_inds1 = filter(x->  !(x in boundary_inds),  sub_nz_inds)
        nz_inds2 = (block_row-1)*block_size .+ sub_nz_inds
        nz_inds3 = (block_row)*block_size .+ sub_nz_inds

        nz_inds = vcat(nz_inds1, nz_inds2, nz_inds3)

    #final row only has two blocks and has boundary conditions.
    elseif block_row == max_block_row

        nz_inds1 = (max_block_row-2)*block_size .+ sub_nz_inds
        nz_inds2 = filter(x->  !(x in boundary_inds), (max_block_row-1)*block_size .+ sub_nz_inds)
        nz_inds = vcat(nz_inds1, nz_inds2)
        
    #second last row also has boundary conditions.
    elseif block_row == max_block_row-1

        nz_inds1 = (max_block_row-3)*block_size .+ sub_nz_inds
        nz_inds2 = (max_block_row-2)*block_size .+ sub_nz_inds
        nz_inds3 = filter(x->  !(x in boundary_inds), (max_block_row-1)*block_size .+ sub_nz_inds)
        nz_inds = vcat(nz_inds1, nz_inds2, nz_inds3)

    #otherwise, there will be 3 full blocks, but there location will depend on block_rows
    else
        nz_inds1 = (block_row-2)*block_size .+ sub_nz_inds
        nz_inds2 = (block_row-1)*block_size .+ sub_nz_inds
        nz_inds3 = (block_row)*block_size .+ sub_nz_inds

        nz_inds = vcat(nz_inds1, nz_inds2, nz_inds3)
    end


    return nz_inds



end


function compute_block_size(grids::MID.FFSGridsT)

    return 4 * grids.θ.N * grids.ζ.count
end

function compute_block_size(grids::MID.FSSGridsT)

    return 2 * grids.θ.count * grids.ζ.count
end

function compute_block_size(grids::MID.FFFGridsT)

    return 8 * grids.θ.N * grids.ζ.N
end

function compute_nz_inds(ind, grids::MID.FFSGridsT, indslocal, boundary_inds)

    #each radial block is made up from Nθ rows of θ blocks.
    #this behaves v similar to FSS case, but we now have periodicity.

    

    block_size = 4 * grids.θ.N * grids.ζ.count
    #block_row = Int64(div(indslocal[ind], block_size, RoundDown)) + 1
    block_row = div(indslocal[ind], block_size) + 1

    sub_block_size = 4 * grids.ζ.count

    max_sub_block_row = grids.θ.N #total sub blocks per row in each block

    max_block_row = grids.r.N

    #sub_block_row = Int64(mod(div(indslocal[ind], sub_block_size, RoundDown), grids.θ.N)) + 1
    sub_block_row = div(rem(indslocal[ind], block_size), sub_block_size) + 1

    #display(Int64(mod(sub_block_row, grids.θ.N)))

    if sub_block_row == 1
        sub_nz_inds1 = 1:2*sub_block_size
        sub_nz_inds2 = (max_sub_block_row-1)*sub_block_size+1:max_sub_block_row * sub_block_size

        sub_nz_inds = vcat(sub_nz_inds1, sub_nz_inds2)

    elseif sub_block_row == max_sub_block_row
        sub_nz_inds1 = 1:sub_block_size
        sub_nz_inds2 = (max_sub_block_row-2)*sub_block_size+1:max_sub_block_row * sub_block_size

        sub_nz_inds = vcat(sub_nz_inds1, sub_nz_inds2)

    else
        #nzs = 3 * sub_block_size

        sub_nz_inds = collect((sub_block_row-2)*sub_block_size+1:(sub_block_row+1)*sub_block_size)
        #nzstart = (sub_block_row-2) * sub_block_size
        #nzend = (sub_block_row+1) * sub_block_size
    end


    if block_row == 1
        nz_inds1 = filter(x->  !(x in boundary_inds),  sub_nz_inds)
        nz_inds2 = block_size .+ sub_nz_inds
        nz_inds = vcat(nz_inds1, nz_inds2)

    elseif block_row == 2
        nz_inds1 = filter(x->  !(x in boundary_inds),  sub_nz_inds)
        nz_inds2 = (block_row-1)*block_size .+ sub_nz_inds
        nz_inds3 = (block_row)*block_size .+ sub_nz_inds

        nz_inds = vcat(nz_inds1, nz_inds2, nz_inds3)

    elseif block_row == max_block_row

        nz_inds1 = (max_block_row-2)*block_size .+ sub_nz_inds
        nz_inds2 = filter(x->  !(x in boundary_inds), (max_block_row-1)*block_size .+ sub_nz_inds)
        nz_inds = vcat(nz_inds1, nz_inds2)
    
    elseif block_row == max_block_row-1

        nz_inds1 = (max_block_row-3)*block_size .+ sub_nz_inds
        nz_inds2 = (max_block_row-2)*block_size .+ sub_nz_inds
        nz_inds3 = filter(x->  !(x in boundary_inds), (max_block_row-1)*block_size .+ sub_nz_inds)
        nz_inds = vcat(nz_inds1, nz_inds2, nz_inds3)

    else
        nz_inds1 = (block_row-2)*block_size .+ sub_nz_inds
        nz_inds2 = (block_row-1)*block_size .+ sub_nz_inds
        nz_inds3 = (block_row)*block_size .+ sub_nz_inds

        nz_inds = vcat(nz_inds1, nz_inds2, nz_inds3)
    end


    return nz_inds



end


#think this is good now as well.
#there is a lot of repeated calculation in here, not sure it matters tho!
#actually allocating fine, test case only has Nζ=2 so it over allocates.
#fine for normal cases though!
function compute_nz_inds(ind, grids::MID.FFFGridsT, indslocal, boundary_inds)

    #each radial block is made up from Nθ rows of θ blocks.
    #this behaves v similar to FSS case, but we now have periodicity.

    
    #largest blocks, each block is for a single radial point.   
    block_size = compute_block_size(grids)

    block_row = div(indslocal[ind], block_size) + 1


    #within each block there are smaller subblocks,
    #each corresponding to a single θ point.
    sub_block_size = 8 * grids.ζ.N

    #sub_sub_block, which is always size 8.
    #finally within the sub blocks there are 8^2 sub-sub-blocks, each corresponding to a single ζ point.
    ss_block_size = 8

    max_ss_block_row = grids.ζ.N


    max_sub_block_row = grids.θ.N #total sub blocks per row in each block

    max_block_row = grids.r.N


    sub_block_row = div(rem(indslocal[ind], block_size), sub_block_size) + 1


    ss_block_row = div(rem(rem(indslocal[ind], block_size), sub_block_size), ss_block_size) + 1

    """
    How this works:
    We find the indicies for the smallest blocks first.
    Then we construct the next size up blocks but combining multiple of the smaller blocks.

    For the two smallest blocks, we have to consider periodicity, meaning for row 1 and row N,
    there will be two blocks at the start (end) and a single block at the end (start).

    For the largest block, for r, we don't need periodicity, so at row 1 and N there will just be two blocks at the start (end). However for this case we need to eliminate the boundary conditions, 
    which is what filter does. Note this is only applied to the appropriate blocks for speed.
    These must be applied at row 1, 2, N, N-1.

    """
    if ss_block_row == 1
        #print("true")
        ss_nz_inds1 = 1:2*ss_block_size
        ss_nz_inds2 = (max_ss_block_row-1)*ss_block_size+1:max_ss_block_row * ss_block_size
        ss_nz_inds = vcat(ss_nz_inds1, ss_nz_inds2)
    elseif ss_block_row == max_ss_block_row
        ss_nz_inds1 = 1:ss_block_size
        ss_nz_inds2 = (max_ss_block_row-2)*ss_block_size+1:max_ss_block_row * ss_block_size
        ss_nz_inds = vcat(ss_nz_inds1, ss_nz_inds2)
    else
        ss_nz_inds = collect((ss_block_row-2)*ss_block_size+1:(ss_block_row+1)*ss_block_size) 
    end



    if sub_block_row == 1
        sub_nz_inds1 = ss_nz_inds
        sub_nz_inds2 = sub_block_size .+ ss_nz_inds

        sub_nz_inds3 = (max_sub_block_row-1)*sub_block_size .+ ss_nz_inds
        #sub_nz_inds2 = (max_sub_block_row-1)*sub_block_size+1:max_sub_block_row * sub_block_size

        sub_nz_inds = vcat(sub_nz_inds1, sub_nz_inds2, sub_nz_inds3)

    elseif sub_block_row == max_sub_block_row
        sub_nz_inds1 = ss_nz_inds
        sub_nz_inds2 = (max_sub_block_row-2)*sub_block_size .+ ss_nz_inds
        sub_nz_inds3 = (max_sub_block_row-1)*sub_block_size .+ ss_nz_inds

        sub_nz_inds = vcat(sub_nz_inds1, sub_nz_inds2, sub_nz_inds3)

    else
        #nzs = 3 * sub_block_size
        sub_nz_inds1 = (sub_block_row-2)*sub_block_size .+ ss_nz_inds
        sub_nz_inds2 = (sub_block_row-1)*sub_block_size .+ ss_nz_inds
        sub_nz_inds3 = (sub_block_row)*sub_block_size .+ ss_nz_inds

        sub_nz_inds = vcat(sub_nz_inds1, sub_nz_inds2, sub_nz_inds3)
        #nzstart = (sub_block_row-2) * sub_block_size
        #nzend = (sub_block_row+1) * sub_block_size
    end



    if block_row == 1
        nz_inds1 = filter(x->  !(x in boundary_inds),  sub_nz_inds)
        #nz_inds1 = sub_nz_inds
        nz_inds2 = block_size .+ sub_nz_inds
        nz_inds = vcat(nz_inds1, nz_inds2)

    elseif block_row == 2
        nz_inds1 = filter(x->  !(x in boundary_inds),  sub_nz_inds)
        #nz_inds1 = sub_nz_inds
        nz_inds2 = (block_row-1)*block_size .+ sub_nz_inds
        nz_inds3 = (block_row)*block_size .+ sub_nz_inds

        nz_inds = vcat(nz_inds1, nz_inds2, nz_inds3)

    elseif block_row == max_block_row

        nz_inds1 = (max_block_row-2)*block_size .+ sub_nz_inds
        nz_inds2 = filter(x->  !(x in boundary_inds), (max_block_row-1)*block_size .+ sub_nz_inds)
        #nz_inds2 = (max_block_row-1)*block_size .+ sub_nz_inds
        nz_inds = vcat(nz_inds1, nz_inds2)
    
    elseif block_row == max_block_row-1

        nz_inds1 = (max_block_row-3)*block_size .+ sub_nz_inds
        nz_inds2 = (max_block_row-2)*block_size .+ sub_nz_inds
        nz_inds3 = filter(x->  !(x in boundary_inds), (max_block_row-1)*block_size .+ sub_nz_inds)
        #nz_inds3 = (max_block_row-1)*block_size .+ sub_nz_inds
        nz_inds = vcat(nz_inds1, nz_inds2, nz_inds3)

    else
        nz_inds1 = (block_row-2)*block_size .+ sub_nz_inds
        nz_inds2 = (block_row-1)*block_size .+ sub_nz_inds
        nz_inds3 = (block_row)*block_size .+ sub_nz_inds

        nz_inds = vcat(nz_inds1, nz_inds2, nz_inds3)
    end

    return nz_inds



end