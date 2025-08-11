"""
    par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFFGridsT)

Constructs the W and I matrices in parallel.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFFGridsT)

    #instantiate the grids into arrays
    x1grid, x2grid, x3grid = inst_grids(grids)

    #initialise the two structs.
    met = MetT()
    B = BFieldT()

    #not sure if this should be combined into 1 or something, focus on getting to work first.
    ξx1, wgx1 = gausslegendre(grids.x1.gp) #same as python!
    ξx2, wgx2 = gausslegendre(grids.x2.gp)
    ξx3, wgx3 = gausslegendre(grids.x3.gp)

    #gets the basis 
    S = hermite_basis(ξx1, ξx2, ξx3)

    ts = ones(size(S.H))

    #creates the trial and test function arrays.
    #these store the basis functions for each derivative
    #and finite elements basis 
    Φ = init_basis_function(grids)
    Ψ = init_basis_function(grids)

    #arrays to store the row, column and data of each matrix element
    rows = Array{Int64}(undef, 0) 
    cols = Array{Int64}(undef, 0) 
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)

    #determines the indicies in the matrices corresponding to the dirichlet boundary conditions for r.
    boundary_inds = compute_boundary_inds(grids)

    #generalised eval problem WΦ = ω^2 I Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices I and W.
    I = local_matrix_size(grids)
    W = local_matrix_size(grids)

    #gets range of indicies in the global matrix the proc owns
    indstart, indend = MatGetOwnershipRange(Wmat)

    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()

    #gets the grid points corresponding to the indicies owned by the proc
    #these are the grid points this proc will compute
    grid_points = matrix_to_grid(indstart, indend, grids)

    #main loop
    for (x1ind, x2ind, x3ind) in grid_points

        #r==N case is computed when r==N-1
        if x1ind == grids.x1.N
            continue
        end

        #takes the local ξ arrays to a global arrays around the grid point.
        x1, x2, x3, dx1, dx2, dx3 = local_to_global(x1ind, x2ind, x3ind, ξx1, ξx2, ξx3, x1grid, x2grid, x3grid) 

        #jacobian of the local to global transformation.
        jac = dx1 * dx2 * dx3 / 8 #following thesis!

        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, B, met, prob, x1, x2, x3, tm)

        #adjust the basis functions to the current coordinates/mode numbers considered.
        create_global_basis!(Φ, S, grids.x2.pf, grids.x3.pf, dx1, dx2, dx3, ts)

        #negatives for conjugate, will assume the phase factor is conjugate as well.
        create_global_basis!(Ψ, S, -grids.x2.pf, -grids.x3.pf, dx1, dx2, dx3, ts)

        #loop over the Hermite elements for the trial function
        for trialx1 in 1:4, trialx2 in 1:4, trialx3 in 1:4

            #determines the matrix index for the trial function
            right_ind = grid_to_index(x1ind, x2ind, x3ind, trialx1, trialx2, trialx3, grids)

            for testx1 in 1:4, testx2 in 1:4, testx3 in 1:4
                
                #and for the test function. 
                left_ind = grid_to_index(x1ind, x2ind, x3ind, testx1, testx2, testx3, grids)
                
                #only check for boundaries if this is true
                #no other i's can possibly give boundaries
                if x1ind==1 || x1ind==grids.x1.N-1


                    if left_ind == right_ind && left_ind in boundary_inds

                            #diagonals for boundary conditions are set to 1.
                        push!(rows, left_ind)
                        push!(cols, right_ind)
                        push!(Wdata, 0.5+0.0im)
                        push!(Idata, 0.5+0.0im)
                            
                            #=
                            push!(rows, left_ind)
                            push!(cols, right_ind)
                            push!(Wdata, 0.0+0.0im)
                            push!(Idata, 0.0+0.0im)
                            =#
                            
                            

                    #otherwise the boundaries are set to zero, which for sparse matrices
                    #is the same as leaving blank.
                    elseif left_ind in boundary_inds
                        continue
                    elseif right_ind in boundary_inds
                        continue
                    #otherwise a regular case for these indicies.
                    else
                        #integrate the local contribution to our matrices.
                        Wsum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], W, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                        Isum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], I, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                        push!(rows, left_ind)
                        push!(cols, right_ind)
                        push!(Wdata, Wsum)
                        push!(Idata, Isum)
                    end
                else

                    
                    #=
                    if left_ind == right_ind 
                        push!(rows, left_ind)
                        push!(cols, right_ind)
                        push!(Wdata, 0.0+0.0im)
                        push!(Idata, 0.0+0.0im)
                    end
                    =#
                    #integrate the local contribution to our matrices.
                    Wsum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], W, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                    Isum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], I, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)
                    
                    push!(rows, left_ind)
                    push!(cols, right_ind)
                    push!(Wdata, Wsum)
                    push!(Idata, Isum)
                    
                end
            end
        end
    end
    
    #add the data to the global matrices.
    set_values!(Wmat, rows, cols, Wdata) 
    set_values!(Imat, rows, cols, Idata) 

end



"""
    par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

Constructs the W and I matrices in parallel with qfm surfaces.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

    #creates the surfaces interpolation struct from the surfaces
    #and a temp struct used in the transformation
    surf_itp, sd = create_surf_itp(surfs)

    #instantiate the grids into arrays
    x1grid, x2grid, x3grid = inst_grids(grids)

    #initialise the two structs for each coordinate set
    tor_met = MetT()
    tor_B = BFieldT()
    qfm_met = MetT()
    qfm_B = BFieldT()

    #not sure if this should be combined into 1 or something, focus on getting to work first.
    ξx1, wgx1 = gausslegendre(grids.x1.gp) #same as python!
    ξx2, wgx2 = gausslegendre(grids.x2.gp)
    ξx3, wgx3 = gausslegendre(grids.x3.gp)

    #gets the basis 
    S = hermite_basis(ξx1, ξx2, ξx3)

    ts = ones(size(S.H))

    #creates the trial and test function arrays.
    #these store the basis functions for each derivative
    #and finite elements basis 
    Φ = init_basis_function(grids)
    Ψ = init_basis_function(grids)

    #arrays to store the row, column and data of each matrix element
    rows = Array{Int64}(undef, 0) 
    cols = Array{Int64}(undef, 0) 
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)

    #determines the indicies in the matrices corresponding to the dirichlet boundary conditions for r.
    boundary_inds = compute_boundary_inds(grids)

    #generalised eval problem WΦ = ω^2 I Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices I and W.
    I = local_matrix_size(grids)
    W = local_matrix_size(grids)

    #gets range of indicies in the global matrix the proc owns
    indstart, indend = MatGetOwnershipRange(Wmat)

    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()
    #initialise strict for storing the Coordinate transformation between (r, θ, ζ) and (s, ϑ, φ)
    CT = CoordTransformT()

    #gets the grid points corresponding to the indicies owned by the proc
    #these are the grid points this proc will compute
    grid_points = matrix_to_grid(indstart, indend, grids)

    #main loop
    for (x1ind, x2ind, x3ind) in grid_points

        #r==N case is computed when r==N-1
        if x1ind == grids.x1.N
            continue
        end

        #takes the local ξ arrays to a global arrays around the grid point.
        x1, x2, x3, dx1, dx2, dx3 = local_to_global(x1ind, x2ind, x3ind, ξx1, ξx2, ξx3, x1grid, x2grid, x3grid) 

        #jacobian of the local to global transformation.
        jac = dx1 * dx2 * dx3 / 8 #following thesis!

        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, tor_B, tor_met, qfm_B, qfm_met, prob, x1, x2, x3, tm, surf_itp, CT, sd)

        #adjust the basis functions to the current coordinates/mode numbers considered.
        create_global_basis!(Φ, S, grids.x2.pf, grids.x3.pf, dx1, dx2, dx3, ts)

        #negatives for conjugate, will assume the phase factor is conjugate as well.
        create_global_basis!(Ψ, S, -grids.x2.pf, -grids.x3.pf, dx1, dx2, dx3, ts)

        #loop over the Hermite elements for the trial function
        for trialx1 in 1:4, trialx2 in 1:4, trialx3 in 1:4

            #determines the matrix index for the trial function
            right_ind = grid_to_index(x1ind, x2ind, x3ind, trialx1, trialx2, trialx3, grids)

            for testx1 in 1:4, testx2 in 1:4, testx3 in 1:4
                
                #and for the test function. 
                left_ind = grid_to_index(x1ind, x2ind, x3ind, testx1, testx2, testx3, grids)
                
                #only check for boundaries if this is true
                #no other i's can possibly give boundaries
                if x1ind==1 || x1ind==grids.x1.N-1


                    if left_ind == right_ind && left_ind in boundary_inds

                        #diagonals for boundary conditions are set to 1.
                        push!(rows, left_ind)
                        push!(cols, right_ind)
                        push!(Wdata, 1.0+0.0im)
                        push!(Idata, 1.0+0.0im)
                        

                    #otherwise the boundaries are set to zero, which for sparse matrices
                    #is the same as leaving blank.
                    elseif left_ind in boundary_inds
                        continue
                    elseif right_ind in boundary_inds
                        continue
                    #otherwise a regular case for these indicies.
                    else
                        #integrate the local contribution to our matrices.
                        Wsum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], W, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                        Isum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], I, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                        push!(rows, left_ind)
                        push!(cols, right_ind)
                        push!(Wdata, Wsum)
                        push!(Idata, Isum)
                    end
                else
                    
                    #integrate the local contribution to our matrices.
                    Wsum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], W, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                    Isum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], I, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)
                    
                    push!(rows, left_ind)
                    push!(cols, right_ind)
                    push!(Wdata, Wsum)
                    push!(Idata, Isum)
                    
                end
            end
        end
    end
    
    #add the data to the global matrices.
    set_values!(Wmat, rows, cols, Wdata) 
    set_values!(Imat, rows, cols, Idata) 

end
