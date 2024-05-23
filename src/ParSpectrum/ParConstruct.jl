
"""
Computes the matrices in parallel. Radial finite element grid is split between worker procs.

# Args
- prob::ProblemT Struct containing the functions and parameters that define the problem we are solving
- grids::GridT Grids to solve over.
"""
function par_construct(; prob::MID.ProblemT, grids::MID.GridsT)
    
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm) #rank of each worker
    nprocs = MPI.Comm_size(comm) #total number of workers including root.
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

    #broadcast this array to each worker.
    MPI.Bcast!(splits, root, comm)

    #now each proc takes there appropriate part of the grid
    r_start = splits[rank+1] + 1
    r_end = splits[rank+2]

    #each worker constructs their part of the matrix.
    rows, cols, Wdata, Idata = worker_construct(prob=prob, grids=grids, r_start=r_start, r_end=r_end)

    return rows, cols, Wdata, Idata

end



"""
Constructs a portion of the total matrix. Each worker is given a chunk of the radial finite elements grid based on r_start and r_end. Otherwise this function is essentially identical to MID.construct.

# Args
- prob::ProblemT Struct containing the functions and parameters that define the problem we are solving
- grids::GridT Grids to solve over.
- r_start::Int64 Start of this workers chunk of the radial grid.
- r_end::Int64 End of this workers chunk of the radial grid.
"""
function worker_construct(; prob::MID.ProblemT, grids::MID.GridsT, r_start::Int64, r_end::Int64)

    
    nθ, mlist, θgrid = MID.spectral_grid(grids.pmd)
    nζ, nlist, ζgrid = MID.spectral_grid(grids.tmd)

    #initialise the two structs.
    met = MID.MetT(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), 0.0, zeros(3))
    B = MID.BFieldT(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))

    ξ, wg = gausslegendre(grids.rd.gp) #same as python!

    #gets the basis 
    H, dH, ddH = MID.hermite_basis(ξ)

    #the trial function
    Φ = zeros(ComplexF64, 9, 4, grids.rd.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 9, 4, grids.rd.gp)   



    #For parallel we don't predefine the array size, as each proc will have a different interaction 
    #with the boundaries, changing the sizes.
    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    boundary_inds = MID.compute_boundary_inds(grids.rd.N, grids.pmd.count, grids.tmd.count, collect(mlist))


    I = zeros(ComplexF64, 9, 9, grids.rd.gp, nθ, nζ)
    W = zeros(ComplexF64, 9, 9, grids.rd.gp, nθ, nζ)


    #this gives a warning but seems to work perfectly
    #and incredibly efficiently!
    #this essentially defines a plane to fft I, (and W as they are the same size), which can be exectued
    #via p * I, this is done in place and seems to be mad efficient.
    p = plan_fft!(I, [4, 5])



    #now we loop through the grid
    #only relevant chunk of radial grid is considered.
    for i in r_start:r_end

        r, dr = MID.local_to_global(i, ξ, grids.rd.grid)

        jac = dr/2 #following thesis!


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

                                
                                left_ind = MID.grid_to_index(i, k2, l2, testsf, grids.pmd.count, grids.tmd.count)

                                #only check for boundaries if this is true
                                #no other i's can possibly give boundaries
                                if i==1 || i==grids.rd.N-1


                                    if left_ind == right_ind && left_ind in boundary_inds

                                        push!(rows, left_ind)
                                        push!(cols, right_ind)

                                        push!(Wdata, 1.0 + 0.0im)
                                        push!(Idata, 1.0 + 0.0im)
                                        
                                    
                                    #otherwise the boundaries are set to zero, which for sparse matrices
                                    #is the same as leaving blank.
                                    elseif left_ind in boundary_inds
                                        continue
                                    elseif right_ind in boundary_inds
                                        continue
                                    #otherwise a regular case for these indicies.
                                    else
                                        push!(rows, left_ind)
                                        push!(cols, right_ind)
                                        

                                        Wsum = @views MID.gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], W[:, :, :, mind, nind], wg, jac, grids.rd.gp)


                                        Isum = @views MID.gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], I[:, :, :, mind, nind], wg, jac, grids.rd.gp)


                                        push!(Wdata, Wsum)
                                        push!(Idata, Isum)
                                        
                                    end
                                else
                                    
                                    push!(rows, left_ind)
                                    push!(cols, right_ind)
                                        

                                    Wsum = @views MID.gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], W[:, :, :, mind, nind], wg, jac, grids.rd.gp)


                                    Isum = @views MID.gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], I[:, :, :, mind, nind], wg, jac, grids.rd.gp)

                                    push!(Wdata, Wsum)
                                    push!(Idata, Isum)
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end


    return rows, cols, Wdata, Idata
end