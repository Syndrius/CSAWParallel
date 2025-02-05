
using Printf
indstart = 0
indend = 3336
ind = 2
indslocal = indstart:indend - 1

Nr=100;
Nθ=20;
Nζ=10;
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=Axel_q, geo=geo); 
rgrid = rfem_grid(N=Nr, gp=5)
θgrid = afem_grid(N=Nθ, pf=2, gp=5);
ζgrid = afem_grid(N=Nζ, pf=-2, gp=5);
grids = init_grids(rgrid, θgrid, ζgrid);

bcinds = MID.Indexing.compute_boundary_inds(grids);

println(bcinds[bcinds .< 10])


compute_nz_inds(ind, grids, indslocal, bcinds, indstart, indend)


function compute_nz_inds(ind, grids::MID.FFFGridsT, indslocal, boundary_inds, indstart, indend)

    #each radial block is made up from Nθ rows of θ blocks.
    #this behaves v similar to FSS case, but we now have periodicity.

    
    #largest blocks, each block is for a single radial point.   
    #block_size = compute_block_size(grids)
    block_size = 8 * grids.θ.N * grids.ζ.N
    
    display(indslocal[ind])

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

    @printf("block_row is %d\n", block_row)
    @printf("sub block_row is %d\n", sub_block_row)
    @printf("ss block_row is %d\n", ss_block_row)

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

    if (indslocal[ind]==320)
        display(ss_block_row)
        display(sub_block_row)
        display(block_row)
        display("ss nz")
        display(ss_nz_inds)
        display("s nz")
        display(sub_nz_inds)

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

    dnnz = length(nz_inds[indstart .< nz_inds .<= indend])
    onnz = length(nz_inds[indstart .>= nz_inds])
    onnz += length(nz_inds[nz_inds .> indend])

    return dnnz, onnz



end