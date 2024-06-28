"""
Computes the matrices in parallel. Radial finite element grid is split between worker procs.

# Args
- prob::ProblemT Struct containing the functions and parameters that define the problem we are solving
- grids::GridT Grids to solve over.
"""
#seems to be working.
function par_construct(W, I; prob::MID.ProblemT, grids::MID.GridsT)
    
    #think this is the same for both types!

    #this is probably the key function where we can offer a more sophisticated pestc split.
    #could probbaly also do with a petsc set value, or perhaps we should set like 4 at a time or something? i.e maybe each i worth of values??
    #think that would be a nice balance between the two approaches.
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm) #rank of each worker
    nprocs = MPI.Comm_size(comm) #total number of workers including root.
    root = 0
    counts = zeros(Int64, nprocs)
    splits = zeros(Int64, nprocs+1)

    #here we split up the radial grid for each proc
    if rank==root
    
        counts_guess = Int64(div(grids.r.N-1, nprocs, RoundDown))
        Remainder    = Int64(grids.r.N-1 - counts_guess*nprocs)
        counts[:]   .= counts_guess
        for i in 1:Remainder
            counts[i] += 1
        end
        splits[1:end-1] = cumsum(append!([0], counts))[1:nprocs] 
        splits[end] = grids.r.N - 1
    end

    #broadcast this array to each worker.
    MPI.Bcast!(splits, root, comm)

    #now each proc takes there appropriate part of the grid
    r_start = splits[rank+1] + 1
    r_end = splits[rank+2]

    #each worker constructs their part of the matrix.
    #rows, cols, Wdata, Idata = old_worker_construct(prob=prob, grids=grids, r_start=r_start, r_end=r_end)
    worker_construct(W, I, prob, grids, r_start, r_end)

    #return rows, cols, Wdata, Idata

end



"""
Constructs a portion of the total matrix. Each worker is given a chunk of the radial finite elements grid based on r_start and r_end. Otherwise this function is essentially identical to MID.construct.

# Args
- prob::ProblemT Struct containing the functions and parameters that define the problem we are solving
- grids::GridT Grids to solve over.
- r_start::Int64 Start of this workers chunk of the radial grid.
- r_end::Int64 End of this workers chunk of the radial grid.
"""
#new version that directly adds to the matrices.
#can probably be significantly improved by making each proc add only to the bits it is storing.
function worker_construct(Wf, If, prob::MID.ProblemT, grids::MID.FSSGridsT, r_start::Int64, r_end::Int64)

    #Wf and If are fkn stupid names.


    #nθ, mlist, θgrid = MID.spectral_grid(grids.pmd)
    #nζ, nlist, ζgrid = MID.spectral_grid(grids.tmd)

    rgrid, Nθ, mlist, θgrid, Nζ, nlist, ζgrid = MID.Misc.instantiate_grids(grids)

    #initialise the two structs.
    met = MID.MetT(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), 0.0, zeros(3))
    B = MID.BFieldT(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))

    ξ, wg = gausslegendre(grids.r.gp) #same as python!

    #gets the basis 
    H, dH, ddH = MID.hermite_basis(ξ)

    #the trial function
    Φ = zeros(ComplexF64, 9, 4, grids.r.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 9, 4, grids.r.gp)   

    #rgrid = MID.construct_rgrid(grids)


    #For parallel we don't predefine the array size, as each proc will have a different interaction 
    #with the boundaries, changing the sizes.
    #rows = Array{Int64}(undef, 0) #ints
    #cols = Array{Int64}(undef, 0) #ints
    #Idata = Array{ComplexF64}(undef, 0)
    #Wdata = Array{ComplexF64}(undef, 0)


    boundary_inds = MID.compute_boundary_inds(grids)


    I = zeros(ComplexF64, 9, 9, grids.r.gp, Nθ, Nζ)
    W = zeros(ComplexF64, 9, 9, grids.r.gp, Nθ, Nζ)


    #this gives a warning but seems to work perfectly
    #and incredibly efficiently!
    #this essentially defines a plane to fft I, (and W as they are the same size), which can be exectued
    #via p * I, this is done in place and seems to be mad efficient.
    p = plan_fft!(W, [4, 5])



    #now we loop through the grid
    #only relevant chunk of radial grid is considered.
    for i in r_start:r_end

        r, dr = MID.local_to_global(i, ξ, rgrid)

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
                        mind = mod(k1-k2 + Nθ, Nθ) + 1
                        nind = mod(l1-l2 + Nζ, Nζ) + 1


                        for trialsf in 1:4

                            right_ind = MID.grid_to_index(i, k1, l1, trialsf, grids)

                            for testsf in 1:4

                                
                                left_ind = MID.grid_to_index(i, k2, l2, testsf, grids)

                                #only check for boundaries if this is true
                                #no other i's can possibly give boundaries
                                if i==1 || i==grids.r.N-1


                                    if left_ind == right_ind && left_ind in boundary_inds

                                        #push!(rows, left_ind)
                                        #push!(cols, right_ind)

                                        #not ideal way of doing this, will be much slower but will save memeory hopefully!
                                        set_values!(Wf, [left_ind], [right_ind], [1.0+0.0im]) 
                                        set_values!(If, [left_ind], [right_ind], [1.0+0.0im]) 

                                        #push!(Wdata, 1.0 + 0.0im)
                                        #push!(Idata, 1.0 + 0.0im)
                                        
                                    
                                    #otherwise the boundaries are set to zero, which for sparse matrices
                                    #is the same as leaving blank.
                                    elseif left_ind in boundary_inds
                                        continue
                                    elseif right_ind in boundary_inds
                                        continue
                                    #otherwise a regular case for these indicies.
                                    else
                                        #push!(rows, left_ind)
                                        #push!(cols, right_ind)
                                        

                                        Wsum = @views MID.gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], W[:, :, :, mind, nind], wg, jac, grids.r.gp)


                                        Isum = @views MID.gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], I[:, :, :, mind, nind], wg, jac, grids.r.gp)

                                        set_values!(Wf, [left_ind], [right_ind], [Wsum]) 
                                        set_values!(If, [left_ind], [right_ind], [Isum]) 
                                        #push!(Wdata, Wsum)
                                        #push!(Idata, Isum)
                                        
                                    end
                                else
                                    
                                    #push!(rows, left_ind)
                                    #push!(cols, right_ind)
                                        

                                    Wsum = @views MID.gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], W[:, :, :, mind, nind], wg, jac, grids.r.gp)


                                    Isum = @views MID.gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], I[:, :, :, mind, nind], wg, jac, grids.r.gp)

                                    #push!(Wdata, Wsum)
                                    #push!(Idata, Isum)
                                    set_values!(Wf, [left_ind], [right_ind], [Wsum]) 
                                    set_values!(If, [left_ind], [right_ind], [Isum]) 
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end


    #return rows, cols, Wdata, Idata
end



"""
    construct(prob::ProblemT, grids::FFSGridsT)

Constructs the two matrices using the WeakForm of the SAW governing equation. Uses Finite elements with cubic Hermite polynomials in r and θ and the fourier spectral method in ζ. Returns two sparse matrices.

### Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::FFSGridT - Grids to solve over.
"""
function worker_construct(Wf, If, prob::MID.ProblemT, grids::MID.FFSGridsT, r_start::Int64, r_end)


    #island being nothing is stupid, either need to separate cases with and without island, or construct an empty (A=0) island if it is nothing.

    

    rgrid, θgrid, Nζ, nlist, ζgrid = MID.Misc.instantiate_grids(grids)


    #initialise the two structs.
    met = MID.MetT(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), 0.0, zeros(3))
    B = MID.BFieldT(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))

    #not sure if this should be combined into 1 or something, focus on getting to work first.
    ξr, wgr = gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = gausslegendre(grids.θ.gp)

    #gets the basis 

    S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ = MID.hermite_basis(ξr, ξθ)

    #the trial function
    #4 is the number of Hermite shape functions
    #10 is Φ and all its relevant derivatives.
    #under new method we don't use the zeroth derivative, so these could be replaced with 9 

    #shape of this will be cooked, expect 4-> 16, unsure if we combined rd and θd yet, leave separate for now.
    Φ = zeros(ComplexF64, 9, 16, grids.r.gp, grids.θ.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 9, 16, grids.r.gp, grids.θ.gp)   


    #generalised eval problem WΦ = ω^2 I Φ
    #Imat = zeros(ComplexF64, 2 * rd.grid_size * pmd.count * tmd.count, 2 * rd.grid_size * pmd.count * tmd.count)
    #Wmat = zeros(ComplexF64, 2 * rd.grid_size * pmd.count * tmd.count, 2 * rd.grid_size * pmd.count * tmd.count)

    #probably possible to know the size of this first.
    #can probably determine the maximum sized int we need based on size of matrix.
    #need to determine the size of these bad bois, they are eating up a bit of time.
    #this seems to have made a minimal difference, but is surely better practise right.
    #arr_length = compute_length(grids.rd.N, grids.pmd.count, grids.tmd.count)

    #arr_count = 1 #this may be wrong... gives error for v small matrix...

    #probably won't know the lengths anymore!
    #rows = Array{Int64}(undef, arr_length) #ints
    #cols = Array{Int64}(undef, arr_length) #ints
    #Idata = Array{ComplexF64}(undef, arr_length)
    #Wdata = Array{ComplexF64}(undef, arr_length)

    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #either need a condition in case m=0 or just an error message.
    boundary_inds = MID.compute_boundary_inds(grids)
    #display(size(boundary_inds))
    #display(boundary_inds)

    #these will hopefully be smaller I think!
    I = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, Nζ)
    W = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, Nζ)


    #this gives a warning but seems to work perfectly
    #and incredibly efficiently!
    #this essentially defines a plane to fft I, (and W as they are the same size), which can be exectued
    #via p * I, this is done in place and seems to be mad efficient.
    #only a single ft now!
    p = plan_fft!(W, [5])

    #Isum = 0.0 + 0.0im
    #Wsum = 0.0 + 0.0im


    #now we loop through the grid

    #will we want a clustered θgrid??? probably not???
    #but we can generalise this function as well if we want it to work with θ, even without clustering!
    #rgrid, θgrid = construct_grids_zf(grids)

    for i in r_start:r_end, j in 1:grids.θ.N #go to N for periodicity!


        r, θ, dr, dθ = MID.local_to_global(i, j, ξr, ξθ, rgrid, θgrid) #wot is θgrid? will need to be constructed I think.

        #Hopefully this is correct!
        jac = dr * dθ / 4 #following thesis!


        #I_and_W!(I, W, B, q_profile, met, compute_met, dens, r, θgrid, ζgrid, δ, isl, R0)

        #hopefully this step will be smaller! but we have twice the loop, so everything else will be longer!
        MID.WeakForm.W_and_I!(W, I, met, B, prob, r, θ, ζgrid)
        #W_tor, I_tor = stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)
        #stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)

        #display(W[:, :, 1, 1, 1])


        #uses the fft plan to take the fft of our two matrices.
        p * W
        p * I


        
        #loop over the fourier components of the trial function
        for (l1, n1) in enumerate(nlist)

            #note we haven't implemented pf for r, seems pointless.
            MID.create_local_basis!(Φ, S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ, grids.θ.pf, n1, dr, dθ)

            for (l2, n2) in enumerate(nlist)

                #negatives for conjugate, will assume the phase factor is conjugate as well.
                MID.create_local_basis!(Ψ, S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ, -grids.θ.pf, -n2, dr, dθ)

                #extract the relevant indicies from the ffted matrices.
                #mind = mod(k1-k2 + nθ, nθ) + 1
                nind = mod(l1-l2 + Nζ, Nζ) + 1


                for trialsf in 1:16
                    
                    #may need a θN or something!
                    right_ind = MID.grid_to_index(i, j, l1, trialsf, grids)

                    for testsf in 1:16
                        #display("testsf")
                        #display(testsf)

                        
                        left_ind = MID.grid_to_index(i, j, l2, testsf, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        
                        if i==1 || i==grids.r.N-1


                            if left_ind == right_ind && left_ind in boundary_inds

                                #rows[arr_count] = left_ind
                                #cols[arr_count] = right_ind
                                #Wdata[arr_count] = 1.0 + 0.0im
                                #Idata[arr_count] = 1.0 + 0.0im

                                set_values!(Wf, [left_ind], [right_ind], [1.0+0.0im]) 
                                set_values!(If, [left_ind], [right_ind], [1.0+0.0im])

                                #push!(rows, left_ind)
                                #push!(cols, right_ind)
                                #push!(Wdata, 1.0 + 0.0im)
                                #push!(Idata, 1.0 + 0.0im)
                                
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
                                #push!(rows, left_ind)
                                #push!(cols, right_ind)
                                

                                #TODO -> should be almost the same, 
                                Wsum = @views MID.gauss_integrate(Ψ[:, testsf, :, :], Φ[:, trialsf, :, :], W[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)


                                Isum = @views MID.gauss_integrate(Ψ[:, testsf, :, :], Φ[:, trialsf, :, :], I[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)

                                #push!(Wdata, Wsum)
                                #push!(Idata, Isum)
                                set_values!(Wf, [left_ind], [right_ind], [Wsum]) 
                                set_values!(If, [left_ind], [right_ind], [Isum])
                            end
                        else
                            
                            #rows[arr_count] = left_ind
                            #cols[arr_count] = right_ind
                            #push!(rows, left_ind)
                            #push!(cols, right_ind)
                                

                            Wsum = @views MID.gauss_integrate(Ψ[:, testsf, :, :], Φ[:, trialsf, :, :], W[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)


                            Isum = @views MID.gauss_integrate(Ψ[:, testsf, :, :], Φ[:, trialsf, :, :], I[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)

                            #push!(Wdata, Wsum)
                            #push!(Idata, Isum)
                            set_values!(Wf, [left_ind], [right_ind], [Wsum]) 
                            set_values!(If, [left_ind], [right_ind], [Isum])
                            
                        end
                    end
                end

            end

        end
    end


    #maybe more consisnt for this function to return the rows and data as per parallal case.
    #Wmat = sparse(rows, cols, Wdata)
    #Imat = sparse(rows, cols, Idata)

    #return Wmat, Imat
end

#note that this function could probably work in general, ideally we would have some kind of if ocndition
#on the type of grids passed in, not sure how that will work though!
#why did this happen.
function par_construct_zf(W, I; prob::MID.ProblemT, grids::MID.GridsT)
    
    
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
    #rows, cols, Wdata, Idata = old_worker_construct(prob=prob, grids=grids, r_start=r_start, r_end=r_end)
    worker_construct_zf(W, I, prob=prob, grids=grids, r_start=r_start, r_end=r_end)

    #return rows, cols, Wdata, Idata

end



"""
Computes the matrices in parallel. Radial finite element grid is split between worker procs.

# Args
- prob::ProblemT Struct containing the functions and parameters that define the problem we are solving
- grids::GridT Grids to solve over.
"""
function old_par_construct(; prob::MID.ProblemT, grids::MID.GridsT)
    
    
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
    rows, cols, Wdata, Idata = old_worker_construct(prob=prob, grids=grids, r_start=r_start, r_end=r_end)
    #worker_construct(W, I, prob=prob, grids=grids, r_start=r_start, r_end=r_end)

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
function old_worker_construct(; prob::MID.ProblemT, grids::MID.GridsT, r_start::Int64, r_end::Int64)

    
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

    rgrid = MID.construct_rgrid(grids)


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

        r, dr = MID.local_to_global(i, ξ, rgrid)

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