
"""
    par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FSSGridsT)

Constructs the W and I matrices in parallel for the fss case.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FSSGridsT)

    #instantiate the grids into arrays.
    x1grid, x2grid, x3grid = inst_grids(grids)

    #for spectral method we need the length of the arrays
    Nx2 = length(x2grid)
    Nx3 = length(x3grid)

    #and the list of modes to consider.
    mlist = mode_list(grids.x2)
    nlist = mode_list(grids.x3)

    #initialise the two structs.
    met = MetT()
    B = BFieldT()

    #compute the gaussian qudrature points for finite elements.
    ξ, wg = gausslegendre(grids.x1.gp) 

    #gets the basis 
    S = hermite_basis(ξ)

    #the trial function
    Φ = zeros(ComplexF64, 4, 9, grids.x1.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 4, 9, grids.x1.gp)   

    #arrays to store the row, column and data of each matrix element
    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)

    #determines the indicies in the matrices corresponding to the dirichlet boundary conditions for x1.
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
    for x1ind in grid_points

        if x1ind == grids.x1.N
            continue
        end

        #takes the local ξ array to a global r array around the grid point.
        x1, dx1 = local_to_global(x1ind, ξ, x1grid)

        #jacobian of the local to global transformation.
        jac = dx1/2 

        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, B, met, prob, x1, x2grid, x3grid, tm)

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
                mind = mod(k1-k2 + Nx2, Nx2) + 1
                nind = mod(l1-l2 + Nx3, Nx3) + 1

                #loop over the Hermite elements for the trial function
                for trialx1 in 1:4

                    #determines the matrix index for the trial function
                    right_ind = grid_to_index(x1ind, k1, l1, trialx1, grids)

                    #and for the test function.
                    for testx1 in 1:4

                        #index for the test function
                        left_ind = grid_to_index(x1ind, k2, l2, testx1, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        if x1ind==1 || x1ind==grids.x1.N-1


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
                                
                                #integrate the local contribution to our matrices.
                                Wsum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], W[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                                Isum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], I[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                                
                            end
                        else
                            
                            #integrate the local contribution to our matrices.
                            Wsum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], W[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                            Isum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], I[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                            push!(rows, left_ind)
                            push!(cols, right_ind)
                            push!(Wdata, Wsum)
                            push!(Idata, Isum)
                            
                        end
                    end
                end
            end
        end
    end

    #combine results into global matrix
    set_values!(Wmat, rows, cols, Wdata) 
    set_values!(Imat, rows, cols, Idata) 

end



"""
    par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FSSGridsT, surfs::Array{QFMSurfaceT})

Constructs the W and I matrices in parallel with qfm surfaces using the spectral method in ϑ, φ.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FSSGridsT, surfs::Array{QFMSurfaceT})

    #creates the surfaces interpolation struct from the surfaces
    #and a temp struct used in the transformation
    surf_itp, sd = create_surf_itp(surfs)

    #instantiate the grids into arrays.
    x1grid, x2grid, x3grid = inst_grids(grids)

    #for spectral method we need the length of the arrays
    Nx2 = length(x2grid)
    Nx3 = length(x3grid)

    #and the list of modes to consider.
    mlist = mode_list(grids.x2)
    nlist = mode_list(grids.x3)

    #initialise the two structs for each coordinate set
    tor_met = MetT()
    tor_B = BFieldT()
    qfm_met = MetT()
    qfm_B = BFieldT()

    #compute the gaussian qudrature points for finite elements.
    ξ, wg = gausslegendre(grids.x1.gp) 

    #gets the basis 
    S = hermite_basis(ξ)

    #the trial function
    Φ = zeros(ComplexF64, 4, 9, grids.x1.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 4, 9, grids.x1.gp)   

    #arrays to store the row, column and data of each matrix element
    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)

    #determines the indicies in the matrices corresponding to the dirichlet boundary conditions for x1.
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
    #initialise strict for storing the Coordinate transformation between (r, θ, ζ) and (s, ϑ, φ)
    CT = CoordTransformT()

    #now we loop through the grid
    for x1ind in grid_points

        if x1ind == grids.x1.N
            continue
        end

        #takes the local ξ array to a global r array around the grid point.
        x1, dx1 = local_to_global(x1ind, ξ, x1grid)

        #jacobian of the local to global transformation.
        jac = dx1/2 

        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, tor_B, tor_met, qfm_B, qfm_met, prob, x1, x2grid, x3grid, tm, surf_itp, CT, sd)

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
                mind = mod(k1-k2 + Nx2, Nx2) + 1
                nind = mod(l1-l2 + Nx3, Nx3) + 1

                #loop over the Hermite elements for the trial function
                for trialx1 in 1:4

                    #determines the matrix index for the trial function
                    right_ind = grid_to_index(x1ind, k1, l1, trialx1, grids)

                    #and for the test function.
                    for testx1 in 1:4

                        #index for the test function
                        left_ind = grid_to_index(x1ind, k2, l2, testx1, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        if x1ind==1 || x1ind==grids.x1.N-1


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
                                
                                #integrate the local contribution to our matrices.
                                Wsum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], W[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                                Isum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], I[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                                
                            end
                        else
                            
                            #integrate the local contribution to our matrices.
                            Wsum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], W[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                            Isum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], I[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                            push!(rows, left_ind)
                            push!(cols, right_ind)
                            push!(Wdata, Wsum)
                            push!(Idata, Isum)
                            
                        end
                    end
                end
            end
        end
    end

    #combine results into global matrix
    set_values!(Wmat, rows, cols, Wdata) 
    set_values!(Imat, rows, cols, Idata) 

end
