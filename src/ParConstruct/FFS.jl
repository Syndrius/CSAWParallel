
"""
    par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFSGridsT)

Construtcs the W and I matrices in parallel.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFSGridsT)



    rgrid, θgrid, ζgrid = inst_grids(grids)

    Nζ = length(ζgrid)
    nlist = mode_list(grids.ζ)

    #initialise the two structs.
    met = MetT()
    B = BFieldT()


    ξr, wgr = gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = gausslegendre(grids.θ.gp)

    #gets the basis 
    S = hermite_basis(ξr, ξθ)


    #creates the trial and test function arrays.
    #these store the basis functions for each derivative
    #and finite elements basis 
    Φ = init_basis_function(grids)
    Ψ = init_basis_function(grids)


    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    boundary_inds = compute_boundary_inds(grids)

    #these will hopefully be smaller I think!
    I = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, Nζ)
    W = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, Nζ)


    #fft plan
    p = plan_fft!(W, [5])

    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()

    indstart, indend = MatGetOwnershipRange(Wmat)


    grid_points = matrix_to_grid(indstart, indend, grids)

    for (rind, θind) in grid_points


        if rind == grids.r.N
            #maybe break? should be last case?
            continue
        end


        r, θ, dr, dθ = local_to_global(rind, θind, ξr, ξθ, rgrid, θgrid) 

        jac = dr * dθ / 4 #following thesis!


        W_and_I!(W, I, B, met, prob, r, θ, ζgrid, tm)


        #uses the fft plan to take the fft of our two matrices.
        p * W
        p * I

        
        #loop over the fourier components of the trial function
        for (l1, n1) in enumerate(nlist)

            #note we haven't implemented pf for r, seems pointless.
            create_local_basis!(Φ, S, grids.θ.pf, n1, dr, dθ)

            for (l2, n2) in enumerate(nlist)

                #negatives for conjugate, will assume the phase factor is conjugate as well.
                create_local_basis!(Ψ, S, -grids.θ.pf, -n2, dr, dθ)

                #extract the relevant indicies from the ffted matrices.
                #mind = mod(k1-k2 + nθ, nθ) + 1
                nind = mod(l1-l2 + Nζ, Nζ) + 1


                for trialr in 1:4, trialθ in 1:4
                    
                    #may need a θN or something!
                    right_ind = grid_to_index(rind, θind, l1, trialr, trialθ, grids)

                    for testr in 1:4, testθ in 1:4
                        
                        left_ind = grid_to_index(rind, θind, l2, testr, testθ, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        
                        if rind==1 || rind==grids.r.N-1


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
                                #rows[arr_count] = left_ind
                                #cols[arr_count] = right_ind
                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                
                                Wsum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], W[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)


                                Isum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], I[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)


                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                                #set_values!(Wmat, [left_ind], [right_ind], [Wsum]) 
                                #set_values!(Imat, [left_ind], [right_ind], [Isum])
                            end
                        else
                            
                            #rows[arr_count] = left_ind
                            #cols[arr_count] = right_ind
                            push!(rows, left_ind)
                            push!(cols, right_ind)
                                


                            Wsum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], W[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)


                            Isum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], I[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)

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
    par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFSGridsT, surfs::Array{QFMSurfaceT})

Constructs the W and I matrices in parallel using qfm surfaces.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFSGridsT, surfs::Array{QFMSurfaceT})


   #haven't changed to s, ϑ, φ out of laziness. 

    rgrid, θgrid, ζgrid = inst_grids(grids)

    Nζ = length(ζgrid)
    nlist = mode_list(grids.ζ)

    #initialise the two structs.
    tor_met = MetT()
    tor_B = BFieldT()
    qfm_met = MetT()
    qfm_B = BFieldT()

    surf_itp, sd = create_surf_itp(surfs)

    #not sure if this should be combined into 1 or something, focus on getting to work first.
    ξr, wgr = gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = gausslegendre(grids.θ.gp)

    #gets the basis 

    S = hermite_basis(ξr, ξθ)


    #shape of this will be cooked, expect 4-> 16, unsure if we combined rd and θd yet, leave separate for now.
    Φ = zeros(ComplexF64, 4, 4, 9, grids.r.gp, grids.θ.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 4, 4, 9, grids.r.gp, grids.θ.gp)   


    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #either need a condition in case m=0 or just an error message.
    boundary_inds = compute_boundary_inds(grids)

    #these will hopefully be smaller I think!
    I = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, Nζ)
    W = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, Nζ)


    p = plan_fft!(W, [5])


    tm = TM()
    CT = CoordTsfmT()
    #now we loop through the grid


    indstart, indend = MatGetOwnershipRange(Wmat)


    grid_points = matrix_to_grid(indstart, indend, grids)

    for (rind, θind) in grid_points

        if rind == grids.r.N
            #maybe break? should be last case?
            continue
        end


        r, θ, dr, dθ = local_to_global(rind, θind, ξr, ξθ, rgrid, θgrid) 

        jac = dr * dθ / 4 #following thesis!

        W_and_I!(W, I, tor_B, tor_met, qfm_B, qfm_met,  prob, r, θ, ζgrid, tm, surf_itp, CT, sd)

        #uses the fft plan to take the fft of our two matrices.
        p * W
        p * I
        
        #loop over the fourier components of the trial function
        for (l1, n1) in enumerate(nlist)

            #note we haven't implemented pf for r, seems pointless.
            create_local_basis!(Φ, S, grids.θ.pf, n1, dr, dθ)

            for (l2, n2) in enumerate(nlist)

                #negatives for conjugate, will assume the phase factor is conjugate as well.
                create_local_basis!(Ψ, S, -grids.θ.pf, -n2, dr, dθ)

                #extract the relevant indicies from the ffted matrices.
                #mind = mod(k1-k2 + nθ, nθ) + 1
                nind = mod(l1-l2 + Nζ, Nζ) + 1


                for trialr in 1:4, trialθ in 1:4
                    
                    #may need a θN or something!
                    right_ind = grid_to_index(rind, θind, l1, trialr, trialθ, grids)

                    for testr in 1:4, testθ in 1:4
                        #display("testsf")
                        #display(testsf)

                        
                        left_ind = grid_to_index(rind, θind, l2, testr, testθ, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        
                        if rind==1 || rind==grids.r.N-1


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
                                #rows[arr_count] = left_ind
                                #cols[arr_count] = right_ind
                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                
                                Wsum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], W[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)


                                Isum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], I[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)


                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                                #set_values!(Wmat, [left_ind], [right_ind], [Wsum]) 
                                #set_values!(Imat, [left_ind], [right_ind], [Isum])
                            end
                        else
                            
                            #rows[arr_count] = left_ind
                            #cols[arr_count] = right_ind
                            push!(rows, left_ind)
                            push!(cols, right_ind)
                                


                            Wsum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], W[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)


                            Isum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], I[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)

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
