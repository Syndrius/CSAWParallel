
#note that these do not work! have not been updated to new version of weak form
#TODO
"""
    construct(prob::ProblemT, grids::FFSGridsT)

Constructs the two matrices using the WeakForm of the SAW governing equation. Uses Finite elements with cubic Hermite polynomials in r and θ and the fourier spectral method in ζ. Returns two sparse matrices.

### Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::FFSGridT - Grids to solve over.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFSGridsT)


    #island being nothing is stupid, either need to separate cases with and without island, or construct an empty (A=0) island if it is nothing.

    

    rgrid, θgrid, ζgrid = inst_grids(grids)

    Nζ = length(ζgrid)
    nlist = mode_list(grids.ζ)

    #initialise the two structs.
    met = metT()
    B = BFieldT()

    #not sure if this should be combined into 1 or something, focus on getting to work first.
    ξr, wgr = gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = gausslegendre(grids.θ.gp)

    #gets the basis 

    S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ = hermite_basis(ξr, ξθ)

    #the trial function
    #4 is the number of Hermite shape functions
    #10 is Φ and all its relevant derivatives.
    #under new method we don't use the zeroth derivative, so these could be replaced with 9 

    #shape of this will be cooked, expect 4-> 16, unsure if we combined rd and θd yet, leave separate for now.
    Φ = zeros(ComplexF64, 4, 4, 9, grids.r.gp, grids.θ.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 4, 4, 9, grids.r.gp, grids.θ.gp)   


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
    boundary_inds = compute_boundary_inds(grids)
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


    indstart, indend = MatGetOwnershipRange(Wmat)

    #grid should be split so that ζ_start=1 and ζ_end = grids.ζ.count
    #r_start, r_end, θ_start, θ_end, _, _ = get_local_grids(indstart, indend, grids)

    #for i in r_start:r_end, j in θ_start:θ_end #go to N for periodicity!
    grid_points = matrix_to_grid(indstart+1, indend, grids)
    #for i in r_start:r_end, j in θ_start:θ_end, k in ζ_start:ζ_end #go to N for periodicity!

    #this is the fix. for splitting the grid up.
    #plus 1 for julia indexing
    #minus 1 on indend as it gives inclusive results.
    for (rind, θind) in grid_points

        #rind, θind, _, _ = index_to_grid(i, grids)

        if rind == grids.r.N
            #maybe break? should be last case?
            continue
        end


        r, θ, dr, dθ = local_to_global(rind, θind, ξr, ξθ, rgrid, θgrid) #wot is θgrid? will need to be constructed I think.

        #Hopefully this is correct!
        jac = dr * dθ / 4 #following thesis!


        #I_and_W!(I, W, B, q_profile, met, compute_met, dens, r, θgrid, ζgrid, δ, isl, R0)

        #hopefully this step will be smaller! but we have twice the loop, so everything else will be longer!
        W_and_I!(W, I, met, B, prob, r, θ, ζgrid)
        #W_tor, I_tor = stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)
        #stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)

        #display(W[:, :, 1, 1, 1])


        #uses the fft plan to take the fft of our two matrices.
        p * W
        p * I


        
        #loop over the fourier components of the trial function
        for (l1, n1) in enumerate(nlist)

            #note we haven't implemented pf for r, seems pointless.
            create_local_basis!(Φ, S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ, grids.θ.pf, n1, dr, dθ)

            for (l2, n2) in enumerate(nlist)

                #negatives for conjugate, will assume the phase factor is conjugate as well.
                create_local_basis!(Ψ, S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ, -grids.θ.pf, -n2, dr, dθ)

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

                                #rows[arr_count] = left_ind
                                #cols[arr_count] = right_ind
                                #Wdata[arr_count] = 1.0 + 0.0im
                                #Idata[arr_count] = 1.0 + 0.0im

                                #set_values!(Wmat, [left_ind], [right_ind], [1.0+0.0im]) 
                                #set_values!(Imat, [left_ind], [right_ind], [1.0+0.0im])

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

    #maybe more consisnt for this function to return the rows and data as per parallal case.
    #Wmat = sparse(rows, cols, Wdata)
    #Imat = sparse(rows, cols, Idata)

    #return Wmat, Imat
end




function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFSGridsT, surfs::Array{QFMSurfaceT})


    #island being nothing is stupid, either need to separate cases with and without island, or construct an empty (A=0) island if it is nothing.

   #haven;t changed to s, ϑ, φ out of laziness. 

    rgrid, θgrid, ζgrid = inst_grids(grids)

    Nζ = length(ζgrid)
    nlist = mode_list(grids.ζ)

    #initialise the two structs.
    tor_met = MetT()
    tor_B = BFieldT()
    qfm_met = MetT()
    qfm_B = BFieldT()

    surf_itp = create_surf_itp(surfs)

    #not sure if this should be combined into 1 or something, focus on getting to work first.
    ξr, wgr = gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = gausslegendre(grids.θ.gp)

    #gets the basis 

    S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ = hermite_basis(ξr, ξθ)

    #the trial function
    #4 is the number of Hermite shape functions
    #10 is Φ and all its relevant derivatives.
    #under new method we don't use the zeroth derivative, so these could be replaced with 9 

    #shape of this will be cooked, expect 4-> 16, unsure if we combined rd and θd yet, leave separate for now.
    Φ = zeros(ComplexF64, 4, 4, 9, grids.r.gp, grids.θ.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 4, 4, 9, grids.r.gp, grids.θ.gp)   


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
    boundary_inds = compute_boundary_inds(grids)
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


    tm = TM()
    CT = CoordTsfmT()
    #now we loop through the grid

    #will we want a clustered θgrid??? probably not???
    #but we can generalise this function as well if we want it to work with θ, even without clustering!
    #rgrid, θgrid = construct_grids_zf(grids)


    indstart, indend = MatGetOwnershipRange(Wmat)

    #grid should be split so that ζ_start=1 and ζ_end = grids.ζ.count
    #r_start, r_end, θ_start, θ_end, _, _ = get_local_grids(indstart, indend, grids)

    #for i in r_start:r_end, j in θ_start:θ_end #go to N for periodicity!
    grid_points = matrix_to_grid(indstart+1, indend, grids)
    #for i in r_start:r_end, j in θ_start:θ_end, k in ζ_start:ζ_end #go to N for periodicity!

    #this is the fix. for splitting the grid up.
    #plus 1 for julia indexing
    #minus 1 on indend as it gives inclusive results.
    for (rind, θind) in grid_points

        #rind, θind, _, _ = index_to_grid(i, grids)

        if rind == grids.r.N
            #maybe break? should be last case?
            continue
        end


        r, θ, dr, dθ = local_to_global(rind, θind, ξr, ξθ, rgrid, θgrid) #wot is θgrid? will need to be constructed I think.

        #Hopefully this is correct!
        jac = dr * dθ / 4 #following thesis!


        #I_and_W!(I, W, B, q_profile, met, compute_met, dens, r, θgrid, ζgrid, δ, isl, R0)

        #hopefully this step will be smaller! but we have twice the loop, so everything else will be longer!
        W_and_I!(W, I, met, B, prob, r, θ, ζgrid)
        #W_tor, I_tor = stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)
        #stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)

        #display(W[:, :, 1, 1, 1])


        #uses the fft plan to take the fft of our two matrices.
        p * W
        p * I


        
        #loop over the fourier components of the trial function
        for (l1, n1) in enumerate(nlist)

            #note we haven't implemented pf for r, seems pointless.
            create_local_basis!(Φ, S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ, grids.θ.pf, n1, dr, dθ)

            for (l2, n2) in enumerate(nlist)

                #negatives for conjugate, will assume the phase factor is conjugate as well.
                create_local_basis!(Ψ, S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ, -grids.θ.pf, -n2, dr, dθ)

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

                                #rows[arr_count] = left_ind
                                #cols[arr_count] = right_ind
                                #Wdata[arr_count] = 1.0 + 0.0im
                                #Idata[arr_count] = 1.0 + 0.0im

                                #set_values!(Wmat, [left_ind], [right_ind], [1.0+0.0im]) 
                                #set_values!(Imat, [left_ind], [right_ind], [1.0+0.0im])

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

    #maybe more consisnt for this function to return the rows and data as per parallal case.
    #Wmat = sparse(rows, cols, Wdata)
    #Imat = sparse(rows, cols, Idata)

    #return Wmat, Imat
end
