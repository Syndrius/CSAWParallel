
#TODO

"""
    construct(prob::ProblemT, grids::FFFGridsT)

Constructs the two matrices using the WeakForm of the SAW governing equation. Uses Finite elements with cubic Hermite polynomials in r, θ and ζ. Returns two sparse matrices.

### Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::FFFGridT - Grids to solve over.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFFGridsT)

    #island being nothing is stupid, either need to separate cases with and without island, or construct an empty (A=0) island if it is nothing.

    
    rgrid, θgrid, ζgrid = inst_grids(grids)


    #initialise the two structs.
    met = MetT()
    B = BFieldT()

    #not sure if this should be combined into 1 or something, focus on getting to work first.
    ξr, wgr = gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = gausslegendre(grids.θ.gp)
    ξζ, wgζ = gausslegendre(grids.ζ.gp)

    #gets the basis 
    #S = hermite_basis(ξr, ξθ, ξζ)
    #ideally these would be combined in some way, this is fkn stupid.
    S = hermite_basis(ξr, ξθ, ξζ)

    #the trial function
    #4 is the number of Hermite shape functions
    #10 is Φ and all its relevant derivatives.
    #under new method we don't use the zeroth derivative, so these could be replaced with 9 

    #order of this is extemely important for @views.
    #because we integrate over each basis individually, the basis indicies (the 4's), should go first,
    #then when we use @views, the entire block of memory is combined for more efficient summation.
    #probably need a proper comment description for this structure as it is cooked beyond belief.
    Φ = zeros(ComplexF64, 4, 4, 4, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 4, 4, 4, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)   


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

    #FFF might actually be the case where we should use these.
    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #either need a condition in case m=0 or just an error message.
    boundary_inds = compute_boundary_inds(grids)
    #display(size(boundary_inds))
    #display(boundary_inds)

    #these will hopefully be smaller I think!
    I = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)
    W = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)


    #this gives a warning but seems to work perfectly
    #and incredibly efficiently!
    #this essentially defines a plane to fft I, (and W as they are the same size), which can be exectued
    #via p * I, this is done in place and seems to be mad efficient.
    #only a single ft now!


    #Isum = 0.0 + 0.0im
    #Wsum = 0.0 + 0.0im


    #so this will work, but really needs ghost cells for this to make sense.
    indstart, indend = MatGetOwnershipRange(Wmat)

    #@printf("I am proc %d, I have %d, %d\n", MPI.Comm_rank(MPI.COMM_WORLD), indstart, indend)

    

    tm = TM()

    #we will start byy allowing the overlap and not including ghost cells.
    #now we loop through the grid

    #plus 1 here fixed issues.
    grid_points = matrix_to_grid(indstart+1, indend, grids)

    #@printf("I am proc %d, I have (%d, %d, %d) and (%d, %d, %d)\n", MPI.Comm_rank(MPI.COMM_WORLD), grid_points[1][1], grid_points[1][2], grid_points[1][3], grid_points[end][1], grid_points[end][2], grid_points[end][3])
    #for i in r_start:r_end, j in θ_start:θ_end, k in ζ_start:ζ_end #go to N for periodicity!

    #this is the fix. for splitting the grid up.
    #plus 1 for julia indexing
    #minus 1 on indend as it gives inclusive results.
    for (rind, θind, ζind) in grid_points
    #for i in indstart+1:indend
        #this is probably wildly inefficeint.
        #this will be called 16 times for the same result...
        #this should at least be quick...
        #we are doing 16 times the same calculation...
        #fkn stupid af.
        #how we do this needs to be improved drastically.
        #rind, θind, ζind, h = MID.index_to_grid(i, grids)

        #v poor fix to prevent us redoing the same calculations.
        #this was defs the problem, this has caused a significant speedup.
        #if h != 1
        #    continue
        #end

        if rind == grids.r.N
            continue
        end


        #display((i, j, k))
        #r, θ, ζ, dr, dθ, dζ = MID.local_to_global(i, j, k, ξr, ξθ, ξζ, rgrid, θgrid, ζgrid) #wot is θgrid? will need to be constructed I think.
        r, θ, ζ, dr, dθ, dζ = local_to_global(rind, θind, ζind, ξr, ξθ, ξζ, rgrid, θgrid, ζgrid) #wot is θgrid? will need to be constructed I think.

        #Hopefully this is correct!
        jac = dr * dθ * dζ / 8 #following thesis!


        #I_and_W!(I, W, B, q_profile, met, compute_met, dens, r, θgrid, ζgrid, δ, isl, R0)

        #hopefully this step will be smaller! but we have twice the loop, so everything else will be longer!
        W_and_I!(W, I, met, B, prob, r, θ, ζ, tm)
        #W_tor, I_tor = stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)
        #stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)

        #display(W[:, :, 1, 1, 1])



        


        #note we haven't implemented pf for r, seems pointless.
        #may be a better way to do this in the future as v little is actually changing in each loop, especiallly since θ and ζ will have a normal grid.
        #ie perhaps this could just be done for each r or something. Not sure this will be the slowdown tho.
        #TODO -> change the inputs to a single structure, this is cooked 
        create_local_basis!(Φ, S, grids.θ.pf, grids.ζ.pf, dr, dθ, dζ)
        #create_local_basis!(Φ, S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ, grids.θ.pf, n1, dr, dθ)


        #negatives for conjugate, will assume the phase factor is conjugate as well.
        create_local_basis!(Ψ, S, -grids.θ.pf, -grids.ζ.pf, dr, dθ, dζ)




        for trialr in 1:4, trialθ in 1:4, trialζ in 1:4

            #display("got to here ok!")
            #may need a θN or something!
            right_ind = grid_to_index(rind, θind, ζind, trialr, trialθ, trialζ, grids)

            for testr in 1:4, testθ in 1:4, testζ in 1:4
                #display("testsf")
                #display(testsf)

                
                left_ind = grid_to_index(rind, θind, ζind, testr, testθ, testζ, grids)


                #only check for boundaries if this is true
                #no other i's can possibly give boundaries
                
                if rind==1 || rind==grids.r.N-1


                    if left_ind == right_ind && left_ind in boundary_inds

                        #rows[arr_count] = left_ind
                        #cols[arr_count] = right_ind
                        #Wdata[arr_count] = 1.0 + 0.0im
                        #Idata[arr_count] = 1.0 + 0.0im
                        #this is happening 4 times...
                        #push!(rows, left_ind)
                        #push!(cols, right_ind)
                        #push!(Wdata, 0.25 + 0.0im)
                        #push!(Idata, 0.25 + 0.0im)

                        #0.25 here is a choice.
                        set_values!(Wmat, [left_ind], [right_ind], [0.25+0.0im])
                        set_values!(Imat, [left_ind], [right_ind], [0.25+0.0im])
                        
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
                        #display("Are we getting stuck here?")
                        #rows[arr_count] = left_ind
                        #cols[arr_count] = right_ind
                        #push!(rows, left_ind)
                        #push!(cols, right_ind)
                        

                        #TODO -> should be almost the same, 
                        Wsum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], W, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)


                        #Isum = @views gauss_integrate(Ψ[:, testr, :, testθ, :, testζ, :], Φ[:, trialr, :, trialθ, :, trialζ, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)
                        Isum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)
                        #Wsum = 1
                        #Isum = 1
                        #display("perhaps?")

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

    #maybe more consisnt for this function to return the rows and data as per parallal case.
    #Wmat = sparse(rows, cols, Wdata)
    #Imat = sparse(rows, cols, Idata)


    #return Wmat, Imat
end



#TODO

"""
    construct(prob::ProblemT, grids::FFFGridsT)

Constructs the two matrices using the WeakForm of the SAW governing equation. Uses Finite elements with cubic Hermite polynomials in r, θ and ζ. Returns two sparse matrices.

### Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::FFFGridT - Grids to solve over.
"""
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

    #island being nothing is stupid, either need to separate cases with and without island, or construct an empty (A=0) island if it is nothing.
    
    #note havent changed to (s, ϑ, φ) out of laziness
    
    rgrid, θgrid, ζgrid = inst_grids(grids)


    #initialise the two structs.
    tor_met = MetT()
    tor_B = BFieldT()
    qfm_met = MetT()
    qfm_B = BFieldT()

    surf_itp = create_surf_itp(surfs)

    #not sure if this should be combined into 1 or something, focus on getting to work first.
    ξr, wgr = gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = gausslegendre(grids.θ.gp)
    ξζ, wgζ = gausslegendre(grids.ζ.gp)

    #gets the basis 
    #S = hermite_basis(ξr, ξθ, ξζ)
    #ideally these would be combined in some way, this is fkn stupid.
    S = hermite_basis(ξr, ξθ, ξζ)

    #the trial function
    #4 is the number of Hermite shape functions
    #10 is Φ and all its relevant derivatives.
    #under new method we don't use the zeroth derivative, so these could be replaced with 9 

    #order of this is extemely important for @views.
    #because we integrate over each basis individually, the basis indicies (the 4's), should go first,
    #then when we use @views, the entire block of memory is combined for more efficient summation.
    #probably need a proper comment description for this structure as it is cooked beyond belief.
    Φ = zeros(ComplexF64, 4, 4, 4, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 4, 4, 4, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)   


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

    #FFF might actually be the case where we should use these.
    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #either need a condition in case m=0 or just an error message.
    boundary_inds = compute_boundary_inds(grids)
    #display(size(boundary_inds))
    #display(boundary_inds)

    #these will hopefully be smaller I think!
    I = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)
    W = zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)


    #this gives a warning but seems to work perfectly
    #and incredibly efficiently!
    #this essentially defines a plane to fft I, (and W as they are the same size), which can be exectued
    #via p * I, this is done in place and seems to be mad efficient.
    #only a single ft now!


    #Isum = 0.0 + 0.0im
    #Wsum = 0.0 + 0.0im


    #so this will work, but really needs ghost cells for this to make sense.
    indstart, indend = MatGetOwnershipRange(Wmat)

    #@printf("I am proc %d, I have %d, %d\n", MPI.Comm_rank(MPI.COMM_WORLD), indstart, indend)

    

    tm = TM()
    CT = CoordTsfmT()

    #we will start byy allowing the overlap and not including ghost cells.
    #now we loop through the grid

    #plus 1 here fixed issues.
    grid_points = matrix_to_grid(indstart+1, indend, grids)

    #@printf("I am proc %d, I have (%d, %d, %d) and (%d, %d, %d)\n", MPI.Comm_rank(MPI.COMM_WORLD), grid_points[1][1], grid_points[1][2], grid_points[1][3], grid_points[end][1], grid_points[end][2], grid_points[end][3])
    #for i in r_start:r_end, j in θ_start:θ_end, k in ζ_start:ζ_end #go to N for periodicity!

    #this is the fix. for splitting the grid up.
    #plus 1 for julia indexing
    #minus 1 on indend as it gives inclusive results.
    for (rind, θind, ζind) in grid_points
    #for i in indstart+1:indend
        #this is probably wildly inefficeint.
        #this will be called 16 times for the same result...
        #this should at least be quick...
        #we are doing 16 times the same calculation...
        #fkn stupid af.
        #how we do this needs to be improved drastically.
        #rind, θind, ζind, h = MID.index_to_grid(i, grids)

        #v poor fix to prevent us redoing the same calculations.
        #this was defs the problem, this has caused a significant speedup.
        #if h != 1
        #    continue
        #end

        if rind == grids.r.N
            continue
        end


        #display((i, j, k))
        #r, θ, ζ, dr, dθ, dζ = MID.local_to_global(i, j, k, ξr, ξθ, ξζ, rgrid, θgrid, ζgrid) #wot is θgrid? will need to be constructed I think.
        r, θ, ζ, dr, dθ, dζ = local_to_global(rind, θind, ζind, ξr, ξθ, ξζ, rgrid, θgrid, ζgrid) #wot is θgrid? will need to be constructed I think.

        #Hopefully this is correct!
        jac = dr * dθ * dζ / 8 #following thesis!


        #I_and_W!(I, W, B, q_profile, met, compute_met, dens, r, θgrid, ζgrid, δ, isl, R0)

        #hopefully this step will be smaller! but we have twice the loop, so everything else will be longer!
        W_and_I!(W, I, tor_met, tor_B, qfm_met, qfm_B, prob, r, θ, ζ, tm, surf_itp, CT)
        #W_tor, I_tor = stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)
        #stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)

        #display(W[:, :, 1, 1, 1])



        


        #note we haven't implemented pf for r, seems pointless.
        #may be a better way to do this in the future as v little is actually changing in each loop, especiallly since θ and ζ will have a normal grid.
        #ie perhaps this could just be done for each r or something. Not sure this will be the slowdown tho.
        #TODO -> change the inputs to a single structure, this is cooked 
        create_local_basis!(Φ, S, grids.θ.pf, grids.ζ.pf, dr, dθ, dζ)
        #create_local_basis!(Φ, S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ, grids.θ.pf, n1, dr, dθ)


        #negatives for conjugate, will assume the phase factor is conjugate as well.
        create_local_basis!(Ψ, S, -grids.θ.pf, -grids.ζ.pf, dr, dθ, dζ)




        for trialr in 1:4, trialθ in 1:4, trialζ in 1:4

            #display("got to here ok!")
            #may need a θN or something!
            right_ind = grid_to_index(rind, θind, ζind, trialr, trialθ, trialζ, grids)

            for testr in 1:4, testθ in 1:4, testζ in 1:4
                #display("testsf")
                #display(testsf)

                
                left_ind = grid_to_index(rind, θind, ζind, testr, testθ, testζ, grids)


                #only check for boundaries if this is true
                #no other i's can possibly give boundaries
                
                if rind==1 || rind==grids.r.N-1


                    if left_ind == right_ind && left_ind in boundary_inds

                        #rows[arr_count] = left_ind
                        #cols[arr_count] = right_ind
                        #Wdata[arr_count] = 1.0 + 0.0im
                        #Idata[arr_count] = 1.0 + 0.0im
                        #this is happening 4 times...
                        #push!(rows, left_ind)
                        #push!(cols, right_ind)
                        #push!(Wdata, 0.25 + 0.0im)
                        #push!(Idata, 0.25 + 0.0im)

                        #0.25 here is a choice.
                        set_values!(Wmat, [left_ind], [right_ind], [0.25+0.0im])
                        set_values!(Imat, [left_ind], [right_ind], [0.25+0.0im])
                        
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
                        #display("Are we getting stuck here?")
                        #rows[arr_count] = left_ind
                        #cols[arr_count] = right_ind
                        #push!(rows, left_ind)
                        #push!(cols, right_ind)
                        

                        #TODO -> should be almost the same, 
                        Wsum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], W, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)


                        #Isum = @views gauss_integrate(Ψ[:, testr, :, testθ, :, testζ, :], Φ[:, trialr, :, trialθ, :, trialζ, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)
                        Isum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)
                        #Wsum = 1
                        #Isum = 1
                        #display("perhaps?")

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

    #maybe more consisnt for this function to return the rows and data as per parallal case.
    #Wmat = sparse(rows, cols, Wdata)
    #Imat = sparse(rows, cols, Idata)


    #return Wmat, Imat
end
