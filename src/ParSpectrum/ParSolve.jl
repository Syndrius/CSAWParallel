"""
Solves generalised eigenvalue problem in parallel using Slepc.

# Args
- rows::Vector{Int64} Array containing the row indexes.
- cols::Vector{Int64} Array containing the column indexes.
- Wdata::Vector{ComplexF64} Array containing the data for W, corresponding to the rows and cols.
- Idata::Vector{ComplexF64} Array containing the data for I, corresponding to the rows and cols.
- σ::Float64 Frequency for 'shift and invert', typically this is the TAE frequency.
- nev::Int64 Number of eigenvalues to solve for.
- n::Int64 Size of the matrix.
- dir::String Directory the results are written to.
"""
function par_solve(W, I, prob::MID.ProblemT, grids::MID.GridsT, dir::String)

    #we will want to postprocess here I think, including normalise etc.
    #no longer think that is a good idea, requires fk loads of extra inputs...
    #so basic functionality is good,
    #but this is completly fked. 
    #Need to:
    # - incorporate into MID
    # - Move to ParPost (probbaly in mid)
    # - Add cases for FFS and FSS 
    # - Make sure this isn't completly cooked for large data (check how much extra time this will add, especially with heaps of cores)
    # - Add phase factor and correct mode labelling.


    
    
    #matrices are then assembled.
    #assemble!(W)
    #assemble!(I)

    #we create the eps object, auto setup uses the slepcargs above.
    eps = create_eps(W, I; auto_setup=true)

    
    #solve the problem.
    #slepc args above will automatically write the solution to file.
    #note sometimes nev is ignored and solution just gives us the number of converged eigenvalues.
    solve!(eps)

    nconv = EPSGetConverged(eps)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        

        #this will store the radial location of the modes.
        #rmax = zeros(Float64, nconv)
        @printf("Number of converged eigenvalues = %d\n", nconv)
        display("Postprocessing...")

        #this will obvs be different for FSS etc.
        #use mkpath to we don't get error if directory already exists
        mkpath(dir*"/efuncs")
        mkpath(dir*"/efuncs_ft")


        #ie eval, rmax, label??
        #label isn't a float though...
        #cont_reconstruction = zeros(Float64, nconv, nconv, nconv)
        mode_loc = zeros(Float64, nconv)
        mlab = zeros(Int64, nconv)
        nlab = zeros(Int64, nconv)
        evals = zeros(ComplexF64, nconv)
    end


    #unclear why this creates two, we only need one.
    #need both of these as hold-over for real scalar_type.
    vecr, veci = MatCreateVecs(W)

    for ieig in 0:nconv-1

        #the extra things returned here are the complex parts,
        #but in scalar_type=complex, the real parts are also comlpex
        #and complex parts are just zero.
        #this is returning the parallel eigenvalue and the parallel eigenvector
        #in Petsc format.
        vp, _, vecp, _ = EPSGetEigenpair(eps, ieig, vecr, veci)

        #display(vpi)
        #this still writes in a stupid formate.
        #viewer = PetscViewerASCIIOpen(MPI.COMM_WORLD, "test.dat")
        #VecView(vecpr, viewer)

        #so this gets a julia array for each proc...
        #create a julia vector for each procs portion of the total vector
        jvecp, _ = VecGetArray(vecp)




        #g_jvec = MPI.Gather(jvec, 20, 0, MPI.COMM_WORLD)
        jvec = MPI.Gather(jvecp, 0, MPI.COMM_WORLD)
        if MPI.Comm_rank(MPI.COMM_WORLD) == 0

            phi, phi_ft, rmode, max_mode = reconstruct(jvec, grids)

            mode_loc[ieig+1] = rmode
            mlab[ieig+1] = max_mode[1]
            nlab[ieig+1] = max_mode[2]

            evals[ieig+1] = vp

            #jdl2 can save julia structs, probably should use that rather than our fkn txt files.
            efunc_label = @sprintf("efun%04d.jld2", ieig+1)


            #efunc_file = dir * "/efuncs/"

            save_object(dir * "/efuncs/"*efunc_label, phi)
            save_object(dir * "/efuncs_ft/"*efunc_label, phi_ft)

            #cont_reconstruction[i]

            #while not ideal, I think this is the way forward, now fourier transform etc and write to file
            #ideally want to get phi structure in file form. 
            #display(jvec)
            #save_object("test.jld2", g_jvec)
        end

        

        #VecViewPETSC_VIEWER_ASCII_COMMON

    end

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0

        #need to normalise...!
        #think this is a terrible structure, but should be pretty small so we can get away with it.
        save(dir*"/cont_reconstruction.jld2", "evals", evals, "rmode", mode_loc, "mlabs", mlab, "nlabs", nlab)
        #save(dir*"/cont_reconstruction.jld2", "rmodes", mode_loc)
        #save(dir*"/cont_reconstruction.jld2", "m labs", mlab)
        #save(dir*"/cont_reconstruction.jld2", "n labs", nlab)
    end


        #open(dir*"")

        #for i in 1:nconv


        #if MPI.Comm_rank(MPI.COMM_WORLD) == 0

            #VecGetArray




    #instead of return eps, want to return a petsc array with all the eigenfunctions inside.
    #-> not sure that will be possible/practical/easy. Cka seems to be creating some cooked data structure.
    #return eps

    


    
    #free the memory used.
    
    destroy!(eps)

    
end


#this should probably be in MID.
#this is probably sub optimal to do this one eval at a time, but we dont often consider many so should be ok.
#this is an awful, suboptimal function, will need to be improved.
function reconstruct(efunc, grids::MID.FFFGridsT)

    phi = zeros(ComplexF64, grids.r.N, grids.θ.N, grids.ζ.N)
    phi_ft = zeros(ComplexF64, grids.r.N, grids.θ.N, grids.ζ.N)


    rgrid, _, _ = MID.instantiate_grids(grids)

    #maxphi = -100
    #we will probably want a plan for fourier transform.


    #think this can be done more efficeintly by knowing the pattern we should be able to set multiple values at once.
    #but this can do for now.
    #eg have a look at cka's way of doing it.
    for i in 1:matrix_dim(grids)

        #note these are the indicies.
        r, θ, ζ, hs = MID.index_to_grid(i, grids)

        #shouldn't need this if we skip over each 8.
        if hs == 1
            #may be the wrong way around!
            #this doesn't seem to have worked as expected tbh!
            phi[r, θ, ζ] = efunc[i] #.* exp(1im * (m * θgrid[θ] + n * ζgrid[ζ]))
        end
    end


    for i in 1:grids.r.N
        #ft in θ and ζ, may want a plan later.
        phi_ft[i, :, :] = fft(phi[i, :, :], [1, 2])
    end 

    rm = zeros(Int64, grids.θ.N, grids.ζ.N)
    ϕm = zeros(Float64, grids.θ.N, grids.ζ.N)

    for j in 1:grids.θ.N, k in 1:grids.ζ.N

        rm[j, k] = argmax(abs.(real.(phi_ft[:, j, k])))
        ϕm[j, k] = abs.(real.(phi_ft[rm[j, k], j, k]))
    end

    max_mode = argmax(ϕm) #still need to convert this to the correct label... 


    #cka does it better than us tbh.
    #max_mode is the index of the maximum mode, essentially useless atm.
    return phi, phi_ft, rgrid[rm[max_mode]], max_mode

end





"""
Solves generalised eigenvalue problem in parallel using Slepc.

# Args
- rows::Vector{Int64} Array containing the row indexes.
- cols::Vector{Int64} Array containing the column indexes.
- Wdata::Vector{ComplexF64} Array containing the data for W, corresponding to the rows and cols.
- Idata::Vector{ComplexF64} Array containing the data for I, corresponding to the rows and cols.
- σ::Float64 Frequency for 'shift and invert', typically this is the TAE frequency.
- nev::Int64 Number of eigenvalues to solve for.
- n::Int64 Size of the matrix.
- dir::String Directory the results are written to.
"""
function old_par_solve(rows::Vector{Int64}, cols::Vector{Int64}, Wdata::Vector{ComplexF64}, Idata::Vector{ComplexF64}; σ::Float64, nev=5::Int64, n::Int64, dir::String)


    
    #eigenvalues are written in matlab format as this offers greater precision than default, 
    #and is easier to read from file
    evals_str = " -eps_view_values :" * dir * "vals.dat:ascii_matlab"
    #eigenfunctions are written as symmodu, which ignores which part of the eigenfunctions is stored by each proc.
    efuncs_str = " -eps_view_vectors :" * dir * "funcs.dat:ascii_symmodu"

    #these are combined with nev and σ, specifiying the setup for slepc.
    #eps_harmonics means we are searching inside the spectrum.
    #-st_type sinvert
    #slepcargs = @sprintf("-eps_nev %d -eps_target %s -eps_harmonic", nev, σ) * evals_str * efuncs_str
    slepcargs = @sprintf("-eps_nev %d -eps_target %s -st_type sinvert -eps_view", nev, σ) * evals_str * efuncs_str
    
    #display(slepcargs)
    SlepcInitialize(slepcargs)

    #The two empty matrices are created.
    W = create_matrix(n, n, auto_setup=true)
    I = create_matrix(n, n, auto_setup=true)

    #Values are added based on Coo format.
    set_values!(W, rows, cols, Wdata)
    set_values!(I, rows, cols, Idata)
    
    #matrices are then assembled.
    assemble!(W)
    assemble!(I)

    #we create the eps object, auto setup uses the slepcargs above.
    eps = create_eps(W, I; auto_setup=true)

    
    #solve the problem.
    #slepc args above will automatically write the solution to file.
    #note sometimes nev is ignored and solution just gives us the number of converged eigenvalues.
    solve!(eps)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display(EPSGetConverged(eps))
    end
    
    #free the memory used.
    destroy!(W)
    destroy!(I)
    destroy!(eps)

    SlepcFinalize()
    
end
