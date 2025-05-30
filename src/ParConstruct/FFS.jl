
"""
    par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFSGridsT)

Construtcs the W and I matrices in parallel.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFSGridsT)

    #instantiate the grids into arrays.
    x1grid, x2grid, x3grid = inst_grids(grids)

    #for spectral method we need the length of the array
    Nx3 = length(x3grid)
    #and the mode list.
    nlist = mode_list(grids.x3)

    #initialise the two structs.
    met = MetT()
    B = BFieldT()

    #compute the gaussian qudrature points for finite elements.
    ξx1, wgx1 = gausslegendre(grids.x1.gp) 
    ξx2, wgx2 = gausslegendre(grids.x2.gp)

    #Gets the Hermite basis for the radial grid and poloidal grid.
    S = hermite_basis(ξx1, ξx2)
    S = hermite_basis(ξx1, ξx2)

    #creates the trial and test function arrays.
    #these store the basis functions for each derivative
    #and finite elements basis 
    Φ = init_basis_function(grids)
    Ψ = init_basis_function(grids)


    #arrays to store the row, column and data of each matrix element
    #used for constructing sparse matrices.
    rows = Array{Int64}(undef, 0) 
    cols = Array{Int64}(undef, 0) 
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #determines the indicies in the matrices corresponding to the dirichlet boundary conditions for x1.
    boundary_inds = compute_boundary_inds(grids)

    #generalised eval problem WΦ = ω^2 I Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices I and W.
    I = local_matrix_size(grids)
    W = local_matrix_size(grids)

    #fft plan
    p = plan_fft!(W, [5])

    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()

    #gets range of indicies in the global matrix the proc owns
    indstart, indend = MatGetOwnershipRange(Wmat)

    #gets the grid points corresponding to the indicies owned by the proc
    #these are the grid points this proc will compute
    grid_points = matrix_to_grid(indstart, indend, grids)


    for (x1ind, x2ind) in grid_points

        if x1ind == grids.x1.N
            #no computation needed at the boundary
            continue
        end

        #takes the local ξ arrays to a global arrays around the grid point.
        x1, x2, dx1, dx2 = local_to_global(x1ind, x2ind, ξx1, ξx2, x1grid, x2grid) 

        #jacobian of the local to global transformation.
        jac = dx1 * dx2 / 4 

        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, B, met, prob, x1, x2, x3grid, tm)

        #uses the fft plan to take the fft of our two matrices.
        p * W
        p * I

        
        #loop over the fourier components of the trial function
        for (l1, n1) in enumerate(nlist)

            #note we haven't implemented pf for x1, seems pointless.
            create_local_basis!(Φ, S, grids.x2.pf, n1, dx1, dx2)

            for (l2, n2) in enumerate(nlist)

                #negatives for conjugate, will assume the phase factor is conjugate as well.
                create_local_basis!(Ψ, S, -grids.x2.pf, -n2, dx1, dx2)

                #extract the relevant indicies from the ffted matrices.
                nind = mod(l1-l2 + Nx3, Nx3) + 1

                #loop over the Hermite elements for the trial function
                for trialx1 in 1:4, trialx2 in 1:4
                    
                    #determines the matrix index for the trial function
                    right_ind = grid_to_index(x1ind, x2ind, l1, trialx1, trialx2, grids)

                    #and for the test function
                    for testx1 in 1:4, testx2 in 1:4
                        
                        #and for the test function
                        left_ind = grid_to_index(x1ind, x2ind, l2, testx1, testx2, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        
                        if x1ind==1 || x1ind==grids.x1.N-1


                            if left_ind == right_ind && left_ind in boundary_inds

                                #diagonals for boundary conditions are set to 1.
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
                                Wsum = @views gauss_integrate(Ψ[testx1, testx2, :, :, :], Φ[trialx1, trialx2, :, :, :], W[:, :, :, :, nind], wgx1, wgx2, jac, grids.x1.gp, grids.x2.gp)

                                Isum = @views gauss_integrate(Ψ[testx1, testx2, :, :, :], Φ[trialx1, trialx2, :, :, :], I[:, :, :, :, nind], wgx1, wgx2, jac, grids.x1.gp, grids.x2.gp)

                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                            end
                        else
                            
                            #integrate the local contribution to our matrices.
                            Wsum = @views gauss_integrate(Ψ[testx1, testx2, :, :, :], Φ[trialx1, trialx2, :, :, :], W[:, :, :, :, nind], wgx1, wgx2, jac, grids.x1.gp, grids.x2.gp)

                            Isum = @views gauss_integrate(Ψ[testx1, testx2, :, :, :], Φ[trialx1, trialx2, :, :, :], I[:, :, :, :, nind], wgx1, wgx2, jac, grids.x1.gp, grids.x2.gp)

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

    #construct the sparse matrix.
    set_values!(Wmat, rows, cols, Wdata) 
    set_values!(Imat, rows, cols, Idata) 

end



"""
    par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFSGridsT, surfs::Array{QFMSurfaceT})

Constructs the W and I matrices in parallel using qfm surfaces.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFSGridsT, surfs::Array{QFMSurfaceT})

    #creates the surfaces interpolation struct from the surfaces
    #and a temp struct used in the transformation
    surf_itp, sd = create_surf_itp(surfs)

    #instantiate the grids into arrays.
    x1grid, x2grid, x3grid = inst_grids(grids)

    #for spectral method we need the length of the array
    Nx3 = length(x3grid)
    #and the mode list.
    nlist = mode_list(grids.x3)

    #initialise the two structs for each coordinate set
    tor_met = MetT()
    tor_B = BFieldT()
    qfm_met = MetT()
    qfm_B = BFieldT()

    #compute the gaussian qudrature points for finite elements.
    ξx1, wgx1 = gausslegendre(grids.x1.gp) 
    ξx2, wgx2 = gausslegendre(grids.x2.gp)

    #Gets the Hermite basis for the radial grid and poloidal grid.
    S = hermite_basis(ξx1, ξx2)
    S = hermite_basis(ξx1, ξx2)

    #creates the trial and test function arrays.
    #these store the basis functions for each derivative
    #and finite elements basis 
    Φ = init_basis_function(grids)
    Ψ = init_basis_function(grids)


    #arrays to store the row, column and data of each matrix element
    #used for constructing sparse matrices.
    rows = Array{Int64}(undef, 0) 
    cols = Array{Int64}(undef, 0) 
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #determines the indicies in the matrices corresponding to the dirichlet boundary conditions for x1.
    boundary_inds = compute_boundary_inds(grids)

    #generalised eval problem WΦ = ω^2 I Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices I and W.
    I = local_matrix_size(grids)
    W = local_matrix_size(grids)

    #fft plan
    p = plan_fft!(W, [5])

    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()
    #initialise strict for storing the Coordinate transformation between (r, θ, ζ) and (s, ϑ, φ)
    CT = CoordTransformT()

    #gets range of indicies in the global matrix the proc owns
    indstart, indend = MatGetOwnershipRange(Wmat)

    #gets the grid points corresponding to the indicies owned by the proc
    #these are the grid points this proc will compute
    grid_points = matrix_to_grid(indstart, indend, grids)


    for (x1ind, x2ind) in grid_points

        if x1ind == grids.x1.N
            #no computation needed at the boundary
            continue
        end

        #takes the local ξ arrays to a global arrays around the grid point.
        x1, x2, dx1, dx2 = local_to_global(x1ind, x2ind, ξx1, ξx2, x1grid, x2grid) 

        #jacobian of the local to global transformation.
        jac = dx1 * dx2 / 4 

        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, tor_B, tor_met, qfm_B, qfm_met,  prob, x1, x2, x3grid, tm, surf_itp, CT, sd)

        #uses the fft plan to take the fft of our two matrices.
        p * W
        p * I

        
        #loop over the fourier components of the trial function
        for (l1, n1) in enumerate(nlist)

            #note we haven't implemented pf for x1, seems pointless.
            create_local_basis!(Φ, S, grids.x2.pf, n1, dx1, dx2)

            for (l2, n2) in enumerate(nlist)

                #negatives for conjugate, will assume the phase factor is conjugate as well.
                create_local_basis!(Ψ, S, -grids.x2.pf, -n2, dx1, dx2)

                #extract the relevant indicies from the ffted matrices.
                nind = mod(l1-l2 + Nx3, Nx3) + 1

                #loop over the Hermite elements for the trial function
                for trialx1 in 1:4, trialx2 in 1:4
                    
                    #determines the matrix index for the trial function
                    right_ind = grid_to_index(x1ind, x2ind, l1, trialx1, trialx2, grids)

                    #and for the test function
                    for testx1 in 1:4, testx2 in 1:4
                        
                        #and for the test function
                        left_ind = grid_to_index(x1ind, x2ind, l2, testx1, testx2, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        
                        if x1ind==1 || x1ind==grids.x1.N-1


                            if left_ind == right_ind && left_ind in boundary_inds

                                #diagonals for boundary conditions are set to 1.
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
                                Wsum = @views gauss_integrate(Ψ[testx1, testx2, :, :, :], Φ[trialx1, trialx2, :, :, :], W[:, :, :, :, nind], wgx1, wgx2, jac, grids.x1.gp, grids.x2.gp)

                                Isum = @views gauss_integrate(Ψ[testx1, testx2, :, :, :], Φ[trialx1, trialx2, :, :, :], I[:, :, :, :, nind], wgx1, wgx2, jac, grids.x1.gp, grids.x2.gp)

                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                            end
                        else
                            
                            #integrate the local contribution to our matrices.
                            Wsum = @views gauss_integrate(Ψ[testx1, testx2, :, :, :], Φ[trialx1, trialx2, :, :, :], W[:, :, :, :, nind], wgx1, wgx2, jac, grids.x1.gp, grids.x2.gp)

                            Isum = @views gauss_integrate(Ψ[testx1, testx2, :, :, :], Φ[trialx1, trialx2, :, :, :], I[:, :, :, :, nind], wgx1, wgx2, jac, grids.x1.gp, grids.x2.gp)

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

    #construct the sparse matrix.
    set_values!(Wmat, rows, cols, Wdata) 
    set_values!(Imat, rows, cols, Idata) 

end
