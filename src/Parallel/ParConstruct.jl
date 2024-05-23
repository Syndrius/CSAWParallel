
#construct the matricies in parallel using MPI.

#not sure if we want file as an input, may be easiest to just always write to the same place then get bash to move it?
function par_construct(; prob::MID.ProblemT, grids::MID.GridsT, file_base="data/")
    
    
    #rd = RDataT(collect(LinRange(0, 1, 50)))
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm) #rank of each worker
    #workers, processors, procs? what the hek
    nprocs = MPI.Comm_size(comm) #total number of workers includeing root.
    #nprocs = 1
    #rank = 0
    root = 0
    counts = zeros(Int64, nprocs)
    splits = zeros(Int64, nprocs+1)

    #here we split up the radial grid for each proc
    if rank==root
    
        counts_guess = Int64(div(grids.rd.N-1, nprocs, RoundDown))
        Remainder    = Int64(grids.rd.N-1 - counts_guess*nprocs)
        counts[:]   .= counts_guess
        for i in 1:Remainder
            counts[i] += 1
        end
        splits[1:end-1] = cumsum(append!([0], counts))[1:nprocs] 
        splits[end] = grids.rd.N - 1
    end

    MPI.Bcast!(splits, root, comm)
    #MPI.Barrier(comm)

    #display(splits)
    #now each proc takes there appropriate part of the grid
    r_start = splits[rank+1] + 1
    r_end = splits[rank+2]
    #println("r_start = $r_start for proc $rank")
    #println("r_end = $r_end for proc $rank")

    rows, cols, Wdata, Idata = worker_construct(prob=prob, grids=grids, r_start=r_start, r_end=r_end)

    #work out how large 
    #counts = Int32.(MPI.Allgather([length(rows)], comm))
    #global_rows = MPI.Gatherv(rows, counts, root, comm)
    #global_cols = MPI.Gatherv(cols, counts, root, comm)
    #global_Wdata = MPI.Gatherv(Wdata, counts, root, comm)
    #global_Idata = MPI.Gatherv(Idata, counts, root, comm)
    
    #maybe, just maybe we can construct the matrix in parallel from these?
    return rows, cols, Wdata, Idata


end


#almost identical to construct, but returns the arrays of data without making sparse, and only works on a fraction of the r grid, as defined by the two indices, r_start and r_end.
function worker_construct(; prob::MID.ProblemT, grids::MID.GridsT, r_start::Int64, r_end::Int64)

    #island being nothing is stupid, either need to separate cases with and without island, or construct an empty (A=0) island if it is nothing.

    nθ, mlist, θgrid = MID.spectral_grid(grids.pmd)
    nζ, nlist, ζgrid = MID.spectral_grid(grids.tmd)


    #initialise the two structs.
    met = MID.MetT(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), 0.0, zeros(3))
    B = MID.BFieldT(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))

    ξ, wg = gausslegendre(grids.rd.gp) #same as python!

    #gets the basis 
    H, dH, ddH = MID.hermite_basis(ξ)

    #the trial function
    #4 is the number of Hermite shape functions
    #10 is Φ and all its relevant derivatives.
    #under new method we don't use the zeroth derivative, so these could be replaced with 9 
    Φ = zeros(ComplexF64, 9, 4, grids.rd.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 9, 4, grids.rd.gp)   


    #generalised eval problem WΦ = ω^2 I Φ
    #Imat = zeros(ComplexF64, 2 * rd.grid_size * pmd.count * tmd.count, 2 * rd.grid_size * pmd.count * tmd.count)
    #Wmat = zeros(ComplexF64, 2 * rd.grid_size * pmd.count * tmd.count, 2 * rd.grid_size * pmd.count * tmd.count)

    #probably possible to know the size of this first.
    #can probably determine the maximum sized int we need based on size of matrix.
    #need to determine the size of these bad bois, they are eating up a bit of time.
    #this seems to have made a minimal difference, but is surely better practise right.
    #this will be wrong again!
    #had a dip, will probably still work if this is to long, but will be a waste of memory!
    #println("N = $(r_end-r_start + 2) for proc ")
    #most procs only need +1, but proc working on right egde needs +2!
    #this is not going to work because different procs will have different changes due to boundary conditions!
    #arr_length = compute_length(r_end-r_start + 2, grids.pmd.count, grids.tmd.count)

    #arr_count = 1

    #rows = Array{Int64}(undef, arr_length) #ints
    #cols = Array{Int64}(undef, arr_length) #ints
    #Idata = Array{ComplexF64}(undef, arr_length)
    #Wdata = Array{ComplexF64}(undef, arr_length)

    #length of the data arrays is a bit all over the shop with mulit procs
    #as some procs will deal with boundaries while other won't
    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #either need a condition in case m=0 or just an error message.
    boundary_inds = MID.compute_boundary_inds(grids.rd.N, grids.pmd.count, grids.tmd.count)
    #display(size(boundary_inds))
    #display(boundary_inds)

    I = zeros(ComplexF64, 9, 9, grids.rd.gp, nθ, nζ)
    W = zeros(ComplexF64, 9, 9, grids.rd.gp, nθ, nζ)


    #this gives a warning but seems to work perfectly
    #and incredibly efficiently!
    #this essentially defines a plane to fft I, (and W as they are the same size), which can be exectued
    #via p * I, this is done in place and seems to be mad efficient.
    p = plan_fft!(I, [4, 5])

    #Isum = 0.0 + 0.0im
    #Wsum = 0.0 + 0.0im

    bounds_count = 0

    #now we loop through the grid

    for i in r_start:r_end

        r, dr = MID.local_to_global(i, ξ, grids.rd.grid)

        jac = dr/2 #following thesis!


        #I_and_W!(I, W, B, q_profile, met, compute_met, dens, r, θgrid, ζgrid, δ, isl, R0)

        MID.WeakForm.W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)


        #uses the fft plan to take the fft of our two matrices.
        p * W
        p * I

        
        #loop over the fourier components of the trial function
        for (k1,m1) in enumerate(mlist)
            for (l1, n1) in enumerate(nlist)

                MID.create_local_basis!(Φ, H, dH, ddH, m1, n1, jac)

                for (k2, m2) in enumerate(mlist)
                
                    for (l2, n2) in enumerate(nlist)

                        #negatives for conjugate
                        MID.create_local_basis!(Ψ, H, dH, ddH, -m2, -n2, jac)

                        #extract the relevant indicies from the ffted matrices.
                        mind = mod(k1-k2 + nθ, nθ) + 1
                        nind = mod(l1-l2 + nζ, nζ) + 1


                        for trialsf in 1:4

                            right_ind = MID.grid_to_index(i, k1, l1, trialsf, grids.pmd.count, grids.tmd.count)

                            for testsf in 1:4
                                #display("testsf")
                                #display(testsf)

                                
                                left_ind = MID.grid_to_index(i, k2, l2, testsf, grids.pmd.count, grids.tmd.count)

                                #only check for boundaries if this is true
                                #no other i's can possibly give boundaries
                                if i==1 || i==grids.rd.N-1


                                    if left_ind == right_ind && left_ind in boundary_inds

                                        #rows[arr_count] = left_ind
                                        #cols[arr_count] = right_ind
                                        #Wdata[arr_count] = 1.0 + 0.0im
                                        #Idata[arr_count] = 1.0 + 0.0im

                                        push!(rows, left_ind)
                                        push!(cols, right_ind)

                                        push!(Wdata, 1.0 + 0.0im)
                                        push!(Idata, 1.0 + 0.0im)
                                        
                                        #arr_count += 1

                                        #bounds_count += 1
                                    
                                    #otherwise the boundaries are set to zero, which for sparse matrices
                                    #is the same as leaving blank.
                                    elseif left_ind in boundary_inds
                                        continue
                                    elseif right_ind in boundary_inds
                                        continue
                                    #otherwise a regular case for these indicies.
                                    else
                                        #rows[arr_count] = left_ind
                                        #cols[arr_count] = right_ind
                                        push!(rows, left_ind)
                                        push!(cols, right_ind)
                                        

                                        Wsum = @views MID.gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], W[:, :, :, mind, nind], wg, jac, grids.rd.gp)


                                        Isum = @views MID.gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], I[:, :, :, mind, nind], wg, jac, grids.rd.gp)

                                        #Wdata[arr_count] = Wsum
                                        #Idata[arr_count] = Isum

                                        push!(Wdata, Wsum)
                                        push!(Idata, Isum)
                                        
                                        #arr_count += 1
                                    end
                                else
                                    
                                    push!(rows, left_ind)
                                    push!(cols, right_ind)
                                    #rows[arr_count] = left_ind
                                    #cols[arr_count] = right_ind
                                        

                                    Wsum = @views MID.gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], W[:, :, :, mind, nind], wg, jac, grids.rd.gp)


                                    Isum = @views MID.gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], I[:, :, :, mind, nind], wg, jac, grids.rd.gp)

                                    #Wdata[arr_count] = Wsum
                                    #Idata[arr_count] = Isum
                                    push!(Wdata, Wsum)
                                    push!(Idata, Isum)
                                    
                                    #arr_count += 1
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    #maybe more consisnt for this function to return the rows and data as per parallal case.
    #Wmat = sparse(rows, cols, Wdata)
    #Imat = sparse(rows, cols, Idata)

    return rows, cols, Wdata, Idata
end