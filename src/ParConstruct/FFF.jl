"""
    par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFFGridsT)

Constructs the W and I matrices in parallel.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFFGridsT)


    
    rgrid, θgrid, ζgrid = inst_grids(grids)


    #initialise the two structs.
    met = MetT()
    B = BFieldT()

    #not sure if this should be combined into 1 or something, focus on getting to work first.
    ξr, wgr = gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = gausslegendre(grids.θ.gp)
    ξζ, wgζ = gausslegendre(grids.ζ.gp)

    #gets the basis 
    S = hermite_basis(ξr, ξθ, ξζ)

    #order of this is extemely important for @views.
    #because we integrate over each basis individually, the basis indicies (the 4's), should go first,
    #then when we use @views, the entire block of memory is combined for more efficient summation.
    #probably need a proper comment description for this structure as it is cooked beyond belief.
    Φ = zeros(ComplexF64, 4, 4, 4, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 4, 4, 4, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)   

    #FFF might actually be the case where we should use these.
    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #either need a condition in case m=0 or just an error message.
    boundary_inds = compute_boundary_inds(grids)

    #these will hopefully be smaller I think!
    I = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)
    W = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)


    indstart, indend = MatGetOwnershipRange(Wmat)

    tm = TM()

    grid_points = matrix_to_grid(indstart, indend, grids)

    for (rind, θind, ζind) in grid_points

        if rind == grids.r.N
            continue
        end


        r, θ, ζ, dr, dθ, dζ = local_to_global(rind, θind, ζind, ξr, ξθ, ξζ, rgrid, θgrid, ζgrid) 

        jac = dr * dθ * dζ / 8 #following thesis!

        W_and_I!(W, I, B, met, prob, r, θ, ζ, tm)


        create_local_basis!(Φ, S, grids.θ.pf, grids.ζ.pf, dr, dθ, dζ)

        #negatives for conjugate, will assume the phase factor is conjugate as well.
        create_local_basis!(Ψ, S, -grids.θ.pf, -grids.ζ.pf, dr, dθ, dζ)


        for trialr in 1:4, trialθ in 1:4, trialζ in 1:4

            right_ind = grid_to_index(rind, θind, ζind, trialr, trialθ, trialζ, grids)

            for testr in 1:4, testθ in 1:4, testζ in 1:4
                
                left_ind = grid_to_index(rind, θind, ζind, testr, testθ, testζ, grids)

                
                if rind==1 || rind==grids.r.N-1


                    if left_ind == right_ind && left_ind in boundary_inds

                        #0.25 here is a choice.
                        set_values!(Wmat, [left_ind], [right_ind], [0.25+0.0im])
                        set_values!(Imat, [left_ind], [right_ind], [0.25+0.0im])
                        

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
                        

                        #TODO -> Can probably make this a bit hekin nicer.
                        Wsum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], W, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)

                        Isum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)

                        set_values!(Wmat, [left_ind], [right_ind], [Wsum])
                        set_values!(Imat, [left_ind], [right_ind], [Isum])

                        #push!(Wdata, Wsum)
                        #push!(Idata, Isum)
                    end
                else
                    
                    #rows[arr_count] = left_ind
                    #cols[arr_count] = right_ind
                    #push!(rows, left_ind)
                    #push!(cols, right_ind)
                        

                    #Wsum = @views gauss_integrate(Ψ[:, testr, :, testθ, :, testζ, :], Φ[:, trialr, :, trialθ, :, trialζ, :], W, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)
                    Wsum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], W, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)


                    Isum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)

                    
                    set_values!(Wmat, [left_ind], [right_ind], [Wsum])
                    set_values!(Imat, [left_ind], [right_ind], [Isum])
                    #Wsum = 1
                    #Isum = 1
                    #push!(Wdata, Wsum)
                    #push!(Idata, Isum)
                    
                end
            end
        end
    end
    

    

    #set_values!(Wmat, rows, cols, Wdata) 
    #set_values!(Imat, rows, cols, Idata) 

end




"""
    par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

Constructs the W and I matrices in parallel with qfm surfaces.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

    
    rgrid, θgrid, ζgrid = inst_grids(grids)

    #initialise the two structs.
    tor_met = MetT()
    tor_B = BFieldT()
    qfm_met = MetT()
    qfm_B = BFieldT()

    surf_itp, sd = create_surf_itp(surfs)

    #not sure if this should be combined into 1 or something, focus on getting to work first.
    ξr, wgr = gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = gausslegendre(grids.θ.gp)
    ξζ, wgζ = gausslegendre(grids.ζ.gp)

    #gets the basis 
    #S = hermite_basis(ξr, ξθ, ξζ)
    #ideally these would be combined in some way, this is fkn stupid.
    S = hermite_basis(ξr, ξθ, ξζ)


    #order of this is extemely important for @views.
    #because we integrate over each basis individually, the basis indicies (the 4's), should go first,
    #then when we use @views, the entire block of memory is combined for more efficient summation.
    #probably need a proper comment description for this structure as it is cooked beyond belief.
    Φ = zeros(ComplexF64, 4, 4, 4, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 4, 4, 4, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)   

    #FFF might actually be the case where we should use these.
    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)

    boundary_inds = compute_boundary_inds(grids)

    #these will hopefully be smaller I think!
    I = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)
    W = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)


    indstart, indend = MatGetOwnershipRange(Wmat)


    tm = TM()
    CT = CoordTsfmT()

    grid_points = matrix_to_grid(indstart, indend, grids)

    for (rind, θind, ζind) in grid_points

        if rind == grids.r.N
            continue
        end


        r, θ, ζ, dr, dθ, dζ = local_to_global(rind, θind, ζind, ξr, ξθ, ξζ, rgrid, θgrid, ζgrid) #wot is θgrid? will need to be constructed I think.

        jac = dr * dθ * dζ / 8 #following thesis!

        W_and_I!(W, I, tor_B, tor_met, qfm_B, qfm_met, prob, r, θ, ζ, tm, surf_itp, CT, sd)

        create_local_basis!(Φ, S, grids.θ.pf, grids.ζ.pf, dr, dθ, dζ)


        #negatives for conjugate, will assume the phase factor is conjugate as well.
        create_local_basis!(Ψ, S, -grids.θ.pf, -grids.ζ.pf, dr, dθ, dζ)

        for trialr in 1:4, trialθ in 1:4, trialζ in 1:4

            right_ind = grid_to_index(rind, θind, ζind, trialr, trialθ, trialζ, grids)

            for testr in 1:4, testθ in 1:4, testζ in 1:4
                
                left_ind = grid_to_index(rind, θind, ζind, testr, testθ, testζ, grids)

                if rind==1 || rind==grids.r.N-1


                    if left_ind == right_ind && left_ind in boundary_inds


                        #0.25 here is a choice.
                        set_values!(Wmat, [left_ind], [right_ind], [0.25+0.0im])
                        set_values!(Imat, [left_ind], [right_ind], [0.25+0.0im])
                    
                    #otherwise the boundaries are set to zero, which for sparse matrices
                    #is the same as leaving blank.
                    elseif left_ind in boundary_inds
                        continue
                    elseif right_ind in boundary_inds
                        continue
                    #otherwise a regular case for these indicies.
                    else

                        #TODO -> should be almost the same, 
                        Wsum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], W, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)


                        Isum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)

                        set_values!(Wmat, [left_ind], [right_ind], [Wsum])
                        set_values!(Imat, [left_ind], [right_ind], [Isum])

                    end
                else
                    
                    Wsum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], W, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)

                    Isum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)

                    
                    set_values!(Wmat, [left_ind], [right_ind], [Wsum])
                    set_values!(Imat, [left_ind], [right_ind], [Isum])
                    
                end
            end
        end
    end
    

    
    #set_values!(Wmat, rows, cols, Wdata) 
    #set_values!(Imat, rows, cols, Idata) 

end
