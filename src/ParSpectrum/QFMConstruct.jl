
#construct for qfm surfaces cases, need to change structure tbh!

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
function par_construct(Wmat::PetscWrap.PetscMat, Imat::PetscWrap.PetscMat, prob::MID.ProblemT, grids::MID.FSSGridsT, surfs::Array{MID.QFMSurfaceT})

    #Wf and If are fkn stupid names.


    #nθ, mlist, θgrid = MID.spectral_grid(grids.pmd)
    #nζ, nlist, ζgrid = MID.spectral_grid(grids.tmd)

    sgrid, ϑgrid, φgrid = inst_grids(grids)

    #inconsistent lables but cbf tbh
    Nθ = length(ϑgrid)
    Nζ = length(φgrid)

    mlist = mode_list(grids.θ)
    nlist = mode_list(grids.ζ)

    #initialise the two structs.
    tor_met = MID.MetT()
    tor_B = MID.BFieldT()
    qfm_met = MID.MetT()
    qfm_B = MID.BFieldT()

    surf_itp = MID.create_surf_itp(surfs)

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


    #this gives a warning but seems to work perfectly
    #and incredibly efficiently!
    #this essentially defines a plane to fft I, (and W as they are the same size), which can be exectued
    #via p * I, this is done in place and seems to be mad efficient.
    p = plan_fft!(W, [4, 5])


    indstart, indend = MatGetOwnershipRange(Wmat)

    #TODO Move this!
    #plus one to shift to julia indexing
    s_start = Int64(indstart / (grids.θ.N * grids.ζ.N * 2)) + 1
    #this has a plus 1 for indexing, and a minus 1 as ownership range returns end+1.
    s_end = Int64(indend / (grids.θ.N * grids.ζ.N * 2))

    @printf("Core %d has %d to %d\n", MPI.Comm_rank(MPI.COMM_WORLD), s_start, s_end)

    if MPI.Comm_rank(MPI.COMM_WORLD) == MPI.Comm_size(MPI.COMM_WORLD)-1
        s_end = s_end-1
    end

    tm = MID.WeakForm.TM()

    CT = MID.CoordTsfmT()
    

    #now we loop through the grid
    #only relevant chunk of radial grid is considered.
    for i in s_start:s_end

        s, ds = local_to_global(i, ξ, sgrid)

        jac = ds/2 #following thesis!


        W_and_I!(W, I, tor_met, tor_B, qfm_met, qfm_B, prob, s, ϑgrid, φgrid, tm, surf_itp, CT)


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

                    right_ind = grid_to_index(i, k1, l1, trialr, grids)

                    for testr in 1:4

                        
                        left_ind = grid_to_index(i, k2, l2, testr, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        if i==1 || i==grids.r.N-1


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
