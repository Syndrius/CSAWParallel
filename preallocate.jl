
#testing some preallocation stuff.
using MID


rank = 0 #rank of each worker
nprocs = 2 #total number of workers including root.
root = 0
counts = zeros(Int64, nprocs)
#splits = zeros(Int64, nprocs+1)

#here we split up the radial grid for each proc
rgrid = init_fem_grid(N=3);
θgrid = init_fem_grid(N=4, pf=2);
#θgrid = init_sm_grid(start=1, count=4)
ζgrid = init_fem_grid(N=4, pf=-2);
#ζgrid = init_sm_grid(start=1, count=3)

grids = init_grids(rgrid, θgrid, ζgrid);

counts_guess = Int64(div(grids.r.N, nprocs, RoundDown))
Remainder    = Int64(grids.r.N - counts_guess*nprocs)
counts[:]   .= counts_guess
for i in 1:Remainder
    counts[i] += 1
end
display(counts)
#splits[1:end-1] = cumsum(append!([0], counts))[1:nprocs] 
#splits[end] = grids.r.N - 1


#MPI.Bcast!(splits, root, comm)


global_n = matrix_dim(grids) 

#if we can do this without knowing the type of grid we can keep this as a single function.
local_n = counts[rank+1] * grids.θ.N * grids.ζ.N * 8





indstart = local_n * rank 
indend = local_n * (rank+1)
#indstart=120
#indend=180
#indend = MatGetOwnershipRange(W)



indslocal = indstart:indend 



dnnz = zeros(Int32, local_n)
onnz = zeros(Int32, local_n)

#block_size = 2 * grids.θ.count * grids.ζ.count

#probbaly going to ignore the boundaries, will be much to complicated,
#and for large grids should become negligible.
#may ignore this for a bit.
boundary_inds = MID.compute_boundary_inds(grids)
boundary_inds = []
#matrix of non-zero inds.
#block_inds = compute_block_inds(grids)



println(dnnz)
println(onnz)

function compute_nz_inds(ind, grids::MID.FFFGridsT)

    #each radial block is made up from Nθ rows of θ blocks.
    #this behaves v similar to FSS case, but we now have periodicity.

    

    block_size = 8 * grids.θ.N * grids.ζ.N
    block_row = Int64(div(indslocal[ind], block_size, RoundDown)) + 1

    #display(block_size)

    sub_block_size = 8 * grids.ζ.N

    #display(sub_block_size)

    #sub_sub_block, which is always size 8.
    ss_block_size = 8

    max_ss_block_row = grids.ζ.N

    #display(max_ss_block_row)


    max_sub_block_row = grids.θ.N #total sub blocks per row in each block

    max_block_row = grids.r.N

    #sub_block_row = Int64(mod(div(indslocal[ind], sub_block_size, RoundDown), grids.θ.N)) + 1
    sub_block_row = mod(mod(indslocal[ind], sub_block_size), grids.θ.N) + 1

    display(sub_block_row)

    #could probably be done in a single mod.
    #ss_block_row = Int64(mod(mod(div(indslocal[ind], ss_block_size, RoundDown), grids.θ.N), grids.ζ.N)) + 1
    #ss_block_row = Int64(mod(div(indslocal[ind], ss_block_size, RoundDown), grids.ζ.N)) + 1

    ss_block_row = mod(mod(indslocal[ind], ss_block_size), grids.ζ.N) + 1

    #display(div(indslocal[ind], ss_block_size))
    #display(Int64(mod(sub_block_row, grids.θ.N)))
    #display(ss_block_row)
    #display(div(indslocal[ind], ss_block_size, RoundDown))
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

    println(ss_nz_inds)


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

    println(sub_nz_inds)


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

    println(nz_inds)
    return nz_inds



end

test_inds = [125, 126] 
for i in test_inds
#for i in 1:local_n

    #may be able to do this all at once at the end??
    #not sure how slow this will be.
    #surely negligible compared with everything else right?
    if indslocal[i]+1 in boundary_inds
        #this only removes the boundaries from the rows not the columns
        #i.e. about half of them.
        #may have to be good enough...
        dnnz[i] = 1
        onnz[i] = 0
        continue
    end 

    nz_inds = compute_nz_inds(i, grids)
        
    dnnz[i] = sum(x-> indstart < x <= indend, nz_inds)
    onnz[i] = sum(x-> x <= indstart, nz_inds)
    onnz[i] += sum(x-> x > indend, nz_inds)
    #onnz[i] = length(nz_inds[indstart .>= nz_inds])
    #onnz[i] += length(nz_inds[nz_inds .> indend]) 

    #if (dnnz[i]==96)
    #    display(i)
    #end

    #dnnz[i] = length(nz_inds[indstart .< nz_inds .<= indend])
    #onnz[i] = length(nz_inds[indstart .>= nz_inds])
    #onnz[i] += length(nz_inds[nz_inds .> indend])

end

println(dnnz)

println(dnnz .+ onnz)


b = [1, 2, 3, 4]
a = [1, 6, 7]

filter(x-> !(x in a), b)
display(a)
display(b)