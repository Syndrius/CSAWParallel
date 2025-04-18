"""
    compute_nz_inds(ind::Int64, grids::FSSGridsT, inds::Array{Int64}, boundary_inds::Array{Int64})

Computes the indicies of the non-zero elements for a given row of the matrix. The non-zeros are split int `blocks` for each radial point. The blocks follow a diagonal pattern. In FSS case, each block is dense.
"""
function compute_nz_inds(ind::Int64, grids::FSSGridsT, inds::UnitRange{Int64}, boundary_inds::Array{Int64})


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



"""
    compute_nz_inds(ind::Int64, grids::FFFGridsT, inds::Array{Int64}, boundary_inds::Array{Int64})

Computes the indicies of the non-zero elements for a given row of the matrix. The non-zeros are split int `blocks` for each radial point. The blocks follow a diagonal pattern. 
"""
function compute_nz_inds(ind::Int64, grids::FFSGridsT, indslocal::UnitRange{Int64}, boundary_inds::Array{Int64})

    #each radial block is made up from Nθ rows of θ blocks.
    #this behaves v similar to FSS case, but we now have periodicity.

    block_size = 4 * grids.θ.N * grids.ζ.N

    block_row = div(indslocal[ind], block_size) + 1

    sub_block_size = 4 * grids.ζ.N

    max_sub_block_row = grids.θ.N #total sub blocks per row in each block

    max_block_row = grids.r.N

    
    sub_block_row = div(rem(indslocal[ind], block_size), sub_block_size) + 1


    if sub_block_row == 1
        sub_nz_inds1 = 1:2*sub_block_size
        sub_nz_inds2 = (max_sub_block_row-1)*sub_block_size+1:max_sub_block_row * sub_block_size

        sub_nz_inds = vcat(sub_nz_inds1, sub_nz_inds2)

    elseif sub_block_row == max_sub_block_row
        sub_nz_inds1 = 1:sub_block_size
        sub_nz_inds2 = (max_sub_block_row-2)*sub_block_size+1:max_sub_block_row * sub_block_size

        sub_nz_inds = vcat(sub_nz_inds1, sub_nz_inds2)

    else

        sub_nz_inds = collect((sub_block_row-2)*sub_block_size+1:(sub_block_row+1)*sub_block_size)
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


"""
    compute_nz_inds(ind::Int64, grids::FFFGridsT, inds::Array{Int64}, boundary_inds::Array{Int64})

Computes the indicies of the non-zero elements for a given row of the matrix. The non-zeros are split int `blocks` for each radial point. The blocks follow a diagonal pattern. 
"""
function compute_nz_inds(ind::Int64, grids::FFFGridsT, indslocal::UnitRange{Int64}, boundary_inds::Array{Int64})

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

        sub_nz_inds = vcat(sub_nz_inds1, sub_nz_inds2, sub_nz_inds3)

    elseif sub_block_row == max_sub_block_row
        sub_nz_inds1 = ss_nz_inds
        sub_nz_inds2 = (max_sub_block_row-2)*sub_block_size .+ ss_nz_inds
        sub_nz_inds3 = (max_sub_block_row-1)*sub_block_size .+ ss_nz_inds

        sub_nz_inds = vcat(sub_nz_inds1, sub_nz_inds2, sub_nz_inds3)
    else
        sub_nz_inds1 = (sub_block_row-2)*sub_block_size .+ ss_nz_inds
        sub_nz_inds2 = (sub_block_row-1)*sub_block_size .+ ss_nz_inds
        sub_nz_inds3 = (sub_block_row)*sub_block_size .+ ss_nz_inds

        sub_nz_inds = vcat(sub_nz_inds1, sub_nz_inds2, sub_nz_inds3)
    end


    if sub_block_row == 1
        sub_nz_inds1 = ss_nz_inds
        sub_nz_inds2 = sub_block_size .+ ss_nz_inds
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


"""
    compute_block_size(grids::FFSGridsT)

Returns the size of each block in the matrix.
"""
function compute_block_size(grids::FFSGridsT)

    return 4 * grids.θ.N * grids.ζ.N
end


"""
    compute_block_size(grids::FSSGridsT)

Returns the size of each block in the matrix.
"""
function compute_block_size(grids::FSSGridsT)

    return 2 * grids.θ.N * grids.ζ.N
end


"""
    compute_block_size(grids::FFFGridsT)

Returns the size of each block in the matrix.
"""
function compute_block_size(grids::FFFGridsT)

    return 8 * grids.θ.N * grids.ζ.N
end
