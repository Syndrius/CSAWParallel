
"""
    par_construct(Pmat::PetscWrap.PetscMat, Qmat::PetscWrap.PetscMat, prob::ProblemT, grids::FSSGridsT)

Constructs the P and Q matrices in parallel for the fss case.
"""
function par_construct(Pmat::PetscWrap.PetscMat, Qmat::PetscWrap.PetscMat, prob::ProblemT, grids::FSSGridsT)

    #instantiate the grids into arrays.
    x1grid, x2grid, x3grid = inst_grids(grids)

    #for spectral method we need the length of the arrays
    Nx2 = length(x2grid)
    Nx3 = length(x3grid)

    #and the list of modes to consider.
    mlist = mode_list(grids.x2)
    nlist = mode_list(grids.x3)

    #initialise the two structs to store the metric and the magnetic field.
    met = MetT()
    B = BFieldT()

    #compute the gaussian quadrature points for finite elements.
    ξ, wg = gauss_points(grids)

    #Gets the Hermite basis for the radial grid.
    S = hermite_basis(ξ)

    #array for storing the scaling of the tangent basis functions when transforming
    #from the local ξ∈[-1, 1] to global x∈[x_i, x_{i+1}] domain
    ts = ones(size(S.H))
    
    #creates the trial and test function arrays.
    #these store the basis functions for each derivative
    #and finite elements basis 
    Φ = init_trial_function(grids)
    Ψ = init_trial_function(grids)
    
    #arrays to store the row, column and data of each matrix element
    #used for constructing sparse matrices.
    rows = Array{Int64}(undef, 0) 
    cols = Array{Int64}(undef, 0) 
    Qdata = Array{ComplexF64}(undef, 0)
    Pdata = Array{ComplexF64}(undef, 0)

    #determines the indicies in the matrices corresponding to the dirichlet boundary conditions for x1.
    boundary_inds = compute_boundary_inds(grids)

    #generalised eval problem PΦ = ω^2 Q Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices Q and P.
    Q = init_local_matrix(grids)
    P = init_local_matrix(grids)

    #creates a fft plan for efficient fft used in spectral method.
    p = plan_fft!(P, [4, 5])

    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()

    #arrays to store the global quadrature points
    x1 = zeros(length(ξ)) 
    Δx = zeros(1)

    #gets the indicies that this core owns
    indstart, indend = MatGetOwnershipRange(Pmat)

    #converts the index range into grid points.
    grid_points = matrix_to_grid(indstart, indend, grids)

    #now we loop through the grid
    for x1ind in grid_points

        if x1ind == grids.x1.N
            continue
        end

        #takes the local ξ array to a global r array around the grid point.
        jac = local_to_global!(x1, Δx, x1ind, ξ, x1grid)

        #computes the contribution to the W and I matrices.
        weak_form!(P, Q, B, met, prob, x1, x2grid, x3grid, tm)

        #uses the fft plan to take the fft of our two matrices.
        p * P
        p * Q
        
        #loop over the fourier components of the trial function
        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)

            #transforms the local basis function to the global.
            update_trial_function!(Φ, S, m1, n1, Δx, ts)

            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)

                #negatives for conjugate in the test function.
                update_trial_function!(Ψ, S, -m2, -n2, Δx, ts)

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

                                push!(Pdata, 1.0 + 0.0im)
                                push!(Qdata, 1.0 + 0.0im)
                            
                            #otherwise the boundaries are set to zero, which for sparse matrices
                            #is the same as leaving blank.
                            elseif left_ind in boundary_inds
                                continue
                            elseif right_ind in boundary_inds
                                continue
                            #otherwise a regular case for these indicies.
                            else
                                
                                #integrate the local contribution to our matrices.
                                Psum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], P[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                                Qsum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], Q[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                push!(Pdata, Psum)
                                push!(Qdata, Qsum)
                                
                            end
                        else
                            
                            #integrate the local contribution to our matrices.
                            Psum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], P[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                            Qsum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], Q[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                            push!(rows, left_ind)
                            push!(cols, right_ind)
                            push!(Pdata, Psum)
                            push!(Qdata, Qsum)
                            
                        end
                    end
                end
            end
        end
    end

    #combine results into global matrix
    set_values!(Pmat, rows, cols, Pdata) 
    set_values!(Qmat, rows, cols, Qdata) 

end



"""
    par_construct(Pmat::PetscWrap.PetscMat, Qmat::PetscWrap.PetscMat, prob::ProblemT, grids::FSSGridsT, surfs::Array{QFMSurfaceT})

Constructs the P and Q matrices in parallel with qfm surfaces using the spectral method in ϑ, φ.
"""
function par_construct(Pmat::PetscWrap.PetscMat, Qmat::PetscWrap.PetscMat, prob::ProblemT, grids::FSSGridsT, surfs::Array{QFMSurfaceT})

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

    #compute the gaussian quadrature points for finite elements.
    ξ, wg = gauss_points(grids)

    #Gets the Hermite basis for the radial grid.
    S = hermite_basis(ξ)

    #array for storing the scaling of the tangent basis functions when transforming
    #from the local ξ∈[-1, 1] to global x∈[x_i, x_{i+1}] domain
    ts = ones(size(S.H))
    
    #creates the trial and test function arrays.
    #these store the basis functions for each derivative
    #and finite elements basis 
    Φ = init_trial_function(grids)
    Ψ = init_trial_function(grids)
    
    #arrays to store the row, column and data of each matrix element
    #used for constructing sparse matrices.
    rows = Array{Int64}(undef, 0) 
    cols = Array{Int64}(undef, 0) 
    Qdata = Array{ComplexF64}(undef, 0)
    Pdata = Array{ComplexF64}(undef, 0)

    #determines the indicies in the matrices corresponding to the dirichlet boundary conditions for x1.
    boundary_inds = compute_boundary_inds(grids)

    #generalised eval problem PΦ = ω^2 Q Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices Q and P.
    Q = init_local_matrix(grids)
    P = init_local_matrix(grids)

    #creates a fft plan for efficient fft used in spectral method.
    p = plan_fft!(P, [4, 5])

    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()

    #initialise strict for storing the Coordinate transformation between (r, θ, ζ) and (s, ϑ, φ)
    CT = CoordTransformT()

    #arrays to store the global quadrature points
    x1 = zeros(length(ξ)) 
    Δx = zeros(1)

    #creates the surfaces interpolation struct from the surfaces
    #and a temp struct used in the transformation
    surf_itp, sd = create_surf_itp(surfs)

    #gets the indicies that this core owns
    indstart, indend = MatGetOwnershipRange(Pmat)

    #converts the index range into grid points.
    grid_points = matrix_to_grid(indstart, indend, grids)

    #now we loop through the grid
    for x1ind in grid_points

        if x1ind == grids.x1.N
            continue
        end

        #takes the local ξ array to a global r array around the grid point.
        jac = local_to_global!(x1, Δx, x1ind, ξ, x1grid)

        #computes the contribution to the W and I matrices.
        weak_form!(P, Q, tor_B, tor_met, qfm_B, qfm_met, prob, x1, x2grid, x3grid, tm, surf_itp, CT, sd)

        #uses the fft plan to take the fft of our two matrices.
        p * P
        p * Q
        
        #loop over the fourier components of the trial function
        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)

            #transforms the local basis function to the global.
            update_trial_function!(Φ, S, m1, n1, Δx, ts)

            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)

                #negatives for conjugate in the test function.
                update_trial_function!(Ψ, S, -m2, -n2, Δx, ts)

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

                                push!(Pdata, 1.0 + 0.0im)
                                push!(Qdata, 1.0 + 0.0im)
                                
                            #otherwise the boundaries are set to zero, which for sparse matrices
                            #is the same as leaving blank.
                            elseif left_ind in boundary_inds
                                continue
                            elseif right_ind in boundary_inds
                                continue
                            #otherwise a regular case for these indicies.
                            else
                                
                                #integrate the local contribution to our matrices.
                                Psum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], P[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                                Qsum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], Q[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                push!(Pdata, Psum)
                                push!(Qdata, Qsum)
                                
                            end
                        else
                            
                            #integrate the local contribution to our matrices.
                            Psum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], P[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                            Qsum = @views gauss_integrate(Ψ[testx1, :, :], Φ[trialx1, :, :], Q[:, :, :, mind, nind], wg, jac, grids.x1.gp)

                            push!(rows, left_ind)
                            push!(cols, right_ind)
                            push!(Pdata, Psum)
                            push!(Qdata, Qsum)
                            
                        end
                    end
                end
            end
        end
    end

    #combine results into global matrix
    set_values!(Pmat, rows, cols, Pdata) 
    set_values!(Qmat, rows, cols, Qdata) 

end
