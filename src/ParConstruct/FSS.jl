
"""
    par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FSSGridsT)

Constructs the W and I matrices in parallel for the fss case.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FSSGridsT)

    #instantiate the grids into arrays.
    rgrid, θgrid, ζgrid = inst_grids(grids)

    #for spectral method we need the length of the arrays
    Nθ = length(θgrid)
    Nζ = length(ζgrid)

    #and the list of modes to consider.
    mlist = mode_list(grids.θ)
    nlist = mode_list(grids.ζ)

    #initialise the two structs.
    met = MetT()
    B = BFieldT()

    #compute the gaussian qudrature points for finite elements.
    ξ, wg = gausslegendre(grids.r.gp) 

    #gets the basis 
    S = hermite_basis(ξ)

    #the trial function
    Φ = zeros(ComplexF64, 4, 9, grids.r.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 4, 9, grids.r.gp)   


    #For parallel we don't predefine the array size, as each proc will have a different interaction 
    #with the boundaries, changing the sizes.
    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)

    boundary_inds = compute_boundary_inds(grids)

    #generalised eval problem WΦ = ω^2 I Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices I and W.
    I = local_matrix_size(grids)
    W = local_matrix_size(grids)


    #plan for the fft
    p = plan_fft!(W, [4, 5])


    #gets the indicies that this core owns
    indstart, indend = MatGetOwnershipRange(Wmat)

    #converts the index range into grid points.
    grid_points = matrix_to_grid(indstart, indend, grids)

    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()


    #now we loop through the grid
    for rind in grid_points

        if rind == grids.r.N
            #maybe break? should be last case?
            continue
        end

        #takes the local ξ array to a global r array around the grid point.
        r, dr = local_to_global(rind, ξ, rgrid)

        #jacobian of the local to global transformation.
        jac = dr/2 

        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, B, met, prob, r, θgrid, ζgrid, tm)

        #uses the fft plan to take the fft of our two matrices.
        p * W
        p * I
        
        #loop over the fourier components of the trial function
        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)

            #adjust the basis functions to the current coordinates/mode numbers considered.
            create_local_basis!(Φ, S, m1, n1, jac)

            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)

                #negatives for conjugate
                create_local_basis!(Ψ, S, -m2, -n2, jac)

                #extract the relevant indicies from the ffted matrices.
                mind = mod(k1-k2 + Nθ, Nθ) + 1
                nind = mod(l1-l2 + Nζ, Nζ) + 1


                for trialr in 1:4

                    #determines the matrix index for the trial function
                    right_ind = grid_to_index(rind, k1, l1, trialr, grids)

                    for testr in 1:4

                        #index for the test function
                        left_ind = grid_to_index(rind, k2, l2, testr, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        if rind==1 || rind==grids.r.N-1


                            if left_ind == right_ind && left_ind in boundary_inds

                                push!(rows, left_ind)
                                push!(cols, right_ind)

                                #not ideal way of doing this, will be much slower but will save memeory hopefully!
                                #set_values!(Wmat, [left_ind], [right_ind], [1.0+0.0im]) 
                                #set_values!(Imat, [left_ind], [right_ind], [1.0+0.0im]) 

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
                                

                                Wsum = @views gauss_integrate(Ψ[testr, :, :], Φ[trialr, :, :], W[:, :, :, mind, nind], wg, jac, grids.r.gp)


                                Isum = @views gauss_integrate(Ψ[:, testr, :], Φ[:, trialr, :], I[:, :, :, mind, nind], wg, jac, grids.r.gp)

                                #set_values!(Wmat, [left_ind], [right_ind], [Wsum]) 
                                #set_values!(Imat, [left_ind], [right_ind], [Isum]) 
                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                                
                            end
                        else
                            
                            push!(rows, left_ind)
                            push!(cols, right_ind)
                                

                            Wsum = @views gauss_integrate(Ψ[testr, :, :], Φ[trialr, :, :], W[:, :, :, mind, nind], wg, jac, grids.r.gp)


                            Isum = @views gauss_integrate(Ψ[testr, :, :], Φ[trialr, :, :], I[:, :, :, mind, nind], wg, jac, grids.r.gp)

                            push!(Wdata, Wsum)
                            push!(Idata, Isum)
                            #set_values!(Wmat, [left_ind], [right_ind], [Wsum]) 
                            #set_values!(Imat, [left_ind], [right_ind], [Isum]) 
                            
                        end
                    end
                end
            end
        end
    end

    set_values!(Wmat, rows, cols, Wdata) 
    set_values!(Imat, rows, cols, Idata) 

end




"""
    par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FSSGridsT, surfs::Array{QFMSurfaceT})

Constructs the W and I matrices in parallel with qfm surfaces using the spectral method in ϑ, φ.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FSSGridsT, surfs::Array{QFMSurfaceT})

    #function copied from non qfm, so inconsistent use of r vs s.


    #instantiate the grids into arrays.
    #note that the inputs are the new coords.
    sgrid, ϑgrid, φgrid = inst_grids(grids)

    #inconsistent lables but cbf tbh
    Nθ = length(ϑgrid)
    Nζ = length(φgrid)

    mlist = mode_list(grids.θ)
    nlist = mode_list(grids.ζ)

    #initialise the two structs.
    tor_met = MetT()
    tor_B = BFieldT()
    qfm_met = MetT()
    qfm_B = BFieldT()

    surf_itp, sd = create_surf_itp(surfs)

    ξ, wg = gausslegendre(grids.r.gp) #same as python!

    #gets the basis 
    S = hermite_basis(ξ)

    #the trial function
    Φ = zeros(ComplexF64, 4, 9, grids.r.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 4, 9, grids.r.gp)   

    #rgrid = MID.construct_rgrid(grids)


    #For parallel we don't predefine the array size, as each proc will have a different interaction 
    #with the boundaries, changing the sizes.
    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    boundary_inds = compute_boundary_inds(grids)


    I = zeros(ComplexF64, 9, 9, grids.r.gp, Nθ, Nζ)
    W = zeros(ComplexF64, 9, 9, grids.r.gp, Nθ, Nζ)


    p = plan_fft!(W, [4, 5])



    #gets the indicies that this core owns
    indstart, indend = MatGetOwnershipRange(Wmat)

    #converts the index range into grid points.
    grid_points = matrix_to_grid(indstart, indend, grids)

    tm = WeakForm.TM()

    CT = CoordTsfmT()
    

    #now we loop through the grid
    #only relevant chunk of radial grid is considered.
    for sind in grid_points

        if sind == grids.r.N
            #maybe break? should be last case?
            continue
        end

        s, ds = local_to_global(sind, ξ, sgrid)

        jac = ds/2 #following thesis!


        W_and_I!(W, I, tor_B, tor_met, qfm_B, qfm_met, prob, s, ϑgrid, φgrid, tm, surf_itp, CT, sd)


        #uses the fft plan to take the fft of our two matrices.
        p * W
        p * I

        
        #loop over the fourier components of the trial function
        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)

            create_local_basis!(Φ, S, m1, n1, jac)

            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)

                #negatives for conjugate
                create_local_basis!(Ψ, S, -m2, -n2, jac)

                #extract the relevant indicies from the ffted matrices.
                mind = mod(k1-k2 + Nθ, Nθ) + 1
                nind = mod(l1-l2 + Nζ, Nζ) + 1


                #should be called trial s tbh!
                for trialr in 1:4

                    right_ind = grid_to_index(sind, k1, l1, trialr, grids)

                    for testr in 1:4

                        
                        left_ind = grid_to_index(sind, k2, l2, testr, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        if sind==1 || sind==grids.r.N-1


                            if left_ind == right_ind && left_ind in boundary_inds

                                push!(rows, left_ind)
                                push!(cols, right_ind)

                                #not ideal way of doing this, will be much slower but will save memeory hopefully!
                                #set_values!(Wmat, [left_ind], [right_ind], [1.0+0.0im]) 
                                #set_values!(Imat, [left_ind], [right_ind], [1.0+0.0im]) 

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
                                

                                Wsum = @views gauss_integrate(Ψ[testr, :, :], Φ[trialr, :, :], W[:, :, :, mind, nind], wg, jac, grids.r.gp)


                                Isum = @views gauss_integrate(Ψ[:, testr, :], Φ[:, trialr, :], I[:, :, :, mind, nind], wg, jac, grids.r.gp)

                                #set_values!(Wmat, [left_ind], [right_ind], [Wsum]) 
                                #set_values!(Imat, [left_ind], [right_ind], [Isum]) 
                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                                
                            end
                        else
                            
                            push!(rows, left_ind)
                            push!(cols, right_ind)
                                

                            Wsum = @views gauss_integrate(Ψ[testr, :, :], Φ[trialr, :, :], W[:, :, :, mind, nind], wg, jac, grids.r.gp)


                            Isum = @views gauss_integrate(Ψ[testr, :, :], Φ[trialr, :, :], I[:, :, :, mind, nind], wg, jac, grids.r.gp)

                            push!(Wdata, Wsum)
                            push!(Idata, Isum)
                            #set_values!(Wmat, [left_ind], [right_ind], [Wsum]) 
                            #set_values!(Imat, [left_ind], [right_ind], [Isum]) 
                            
                        end
                    end
                end
            end
        end
    end

    set_values!(Wmat, rows, cols, Wdata) 
    set_values!(Imat, rows, cols, Idata) 


    #return rows, cols, Wdata, Idata
end
