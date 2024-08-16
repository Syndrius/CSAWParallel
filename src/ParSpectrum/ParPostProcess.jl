
"""
    function par_post_process(eps::SlepcWrap.SlepcEPS, dir::String, grids::MID.GridsT, geo::MID.GeoParamsT, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32)

Processes the eigenvalues and eigenfunctions. This includes normalising the eigenvalues, converting the eigenfunctions back to the 3d grid, computing the continuum reconstruction and writing the outputs to file.
"""
function par_post_process(eps::SlepcWrap.SlepcEPS, dir::String, grids::MID.GridsT, geo::MID.GeoParamsT, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32)
    #unclear if this can be done in parallel or not.
    #looks like it may be possible???
    #think start by just getting root to do it all, then maybe we want to split up the solutions.
    

    #will need different versions of this for different grids.
    #looks like a single core is the way.
    #may be possible to parallise later but may not be worth it.

    #number of solutions should be comparable to work of single core anyway.
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0

        #this will obvs be different for FSS etc.
        #use mkpath to we don't get error if directory already exists
        mkpath(dir*"/efuncs")
        mkpath(dir*"/efuncs_ft")


        #ie eval, rmax, label??
        #label isn't a float though...
        #cont_reconstruction = zeros(Float64, nconv, nconv, nconv)
        rmode = zeros(Float64, nconv)
        #mlab = zeros(Int64, nconv)
        #nlab = zeros(Int64, nconv)
        ω = zeros(ComplexF64, nconv)
        modelabs = Tuple{Int, Int}[]
    end


    #unclear why this creates two, we only need one.
    #need both of these as hold-over for real scalar_type.
    

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

            #this allows this to be single function,
            #different cases are covered within this function.
            phi, phi_ft, rm, modelab = process_efunc(jvec, grids)

            


            push!(modelabs, modelab)

            rmode[ieig+1] = rm
            #mlab[ieig+1] = max_mode[1]
            #nlab[ieig+1] = max_mode[2]

            #normalise the eigenvalues
            ω[ieig+1] = geo.R0 * sqrt(vp)

            #jdl2 can save julia structs, probably should use that rather than our fkn txt files.
            efunc_label = @sprintf("efunc%04d.jld2", ieig+1)


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

        evals = MID.EvalsT(ω, rmode, modelabs)

        save_object(dir*"/evals.jld2", evals)

        #need to normalise...!
        #think this is a terrible structure, but should be pretty small so we can get away with it.
        #save(dir*"/evals.jld2", "evals", evals, "rmode", mode_loc, "mlabs", mlab, "nlabs", nlab)
        #save(dir*"/cont_reconstruction.jld2", "rmodes", mode_loc)
        #save(dir*"/cont_reconstruction.jld2", "m labs", mlab)
        #save(dir*"/cont_reconstruction.jld2", "n labs", nlab)
    end

    
end


function process_efunc(efunc, grids::MID.FSSGridsT)

    #this will just be zeros, need a better solution...
    #probably want to un-fourier transform this.
    phi = zeros(ComplexF64, grids.r.N, grids.θ.count, grids.ζ.count)
    phi_ft = zeros(ComplexF64, grids.r.N, grids.θ.count, grids.ζ.count)


    #doing this for every single efunc is probbaly stupid af.
    rgrid, _, mlist, _, _, nlist, _= instantiate_grids(grids)

    #maxphi = -100
    #we will probably want a plan for fourier transform.


    #think this can be done more efficeintly by knowing the pattern we should be able to set multiple values at once.
    #but this can do for now.
    #eg have a look at cka's way of doing it.
    for i in 1:2:matrix_dim(grids)

        #note these are the indicies.
        r, θ, ζ, hs = MID.index_to_grid(i, grids)

        #shouldn't need this if we skip over each 8.
        if hs == 1
            #this is actually already the ft version.
            phi_ft[r, θ, ζ] = efunc[i] #.* exp(1im * (m * θgrid[θ] + n * ζgrid[ζ]))
        end
    end


    #for i in 1:grids.r.N
        #ft in θ and ζ, may want a plan later.
    #    phi_ft[i, :, :] = fft(phi[i, :, :], [1, 2])
    #end 

    rm = zeros(Int64, grids.θ.count, grids.ζ.count)
    ϕm = zeros(Float64, grids.θ.count, grids.ζ.count)

    for j in 1:grids.θ.count, k in 1:grids.ζ.count

        rm[j, k] = argmax(abs.(real.(phi_ft[:, j, k])))
        ϕm[j, k] = abs.(real.(phi_ft[rm[j, k], j, k]))
    end

    max_mode = argmax(ϕm) #still need to convert this to the correct label... 

    mlab = mlist[max_mode[1]]
    nlab = nlist[max_mode[2]]
    


    #cka does it better than us tbh.
    #max_mode is the index of the maximum mode, essentially useless atm.
    return phi, phi_ft, rgrid[rm[max_mode]], (mlab, nlab)

end


function process_efunc(efunc, grids::MID.FFSGridsT)

    phi = zeros(ComplexF64, grids.r.N, grids.θ.N, grids.ζ.count)
    phi_ft = zeros(ComplexF64, grids.r.N, grids.θ.N, grids.ζ.count)


    rgrid, _, _, nlist, _= instantiate_grids(grids)

    #maxphi = -100
    #we will probably want a plan for fourier transform.


    #think this can be done more efficeintly by knowing the pattern we should be able to set multiple values at once.
    #but this can do for now.
    #eg have a look at cka's way of doing it.
    for i in 1:8:matrix_dim(grids)

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
        for n in 1:grids.ζ.count 
            #ft in θ and ζ, may want a plan later.
            phi_ft[i, :, n] = fft(phi[i, :, n], [1])
        end
    end 

    rm = zeros(Int64, grids.θ.N, grids.ζ.count)
    ϕm = zeros(Float64, grids.θ.N, grids.ζ.count)

    for j in 1:grids.θ.N, k in 1:grids.ζ.count

        rm[j, k] = argmax(abs.(real.(phi_ft[:, j, k])))
        ϕm[j, k] = abs.(real.(phi_ft[rm[j, k], j, k]))
    end

    max_mode = argmax(ϕm) #still need to convert this to the correct label... 
    mlab = MID.mode_label(max_mode[1], grids.θ)
    nlab = nlist[max_mode[2]]

    #cka does it better than us tbh.
    #max_mode is the index of the maximum mode, essentially useless atm.
    return phi, phi_ft, rgrid[rm[max_mode]], (mlab, nlab)

end



#this should probably be in MID.
#this is probably sub optimal to do this one eval at a time, but we dont often consider many so should be ok.
#this is an awful, suboptimal function, will need to be improved.
#same as reconstruct phi in MID, but only does a single efunc at a time.
#can probably combine and put into MID. Will never be used by MID though
function process_efunc(efunc, grids::MID.FFFGridsT)

    phi = zeros(ComplexF64, grids.r.N, grids.θ.N, grids.ζ.N)
    phi_ft = zeros(ComplexF64, grids.r.N, grids.θ.N, grids.ζ.N)


    rgrid, _, _ = MID.instantiate_grids(grids)

    #maxphi = -100
    #we will probably want a plan for fourier transform.


    #think this can be done more efficeintly by knowing the pattern we should be able to set multiple values at once.
    #but this can do for now.
    #eg have a look at cka's way of doing it.
    for i in 1:8:matrix_dim(grids)

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

    mlab = mod(max_mode[1]-1, grids.θ.N)
    nlab = mod(max_mode[2]-1, grids.ζ.N)
    if mlab > grids.θ.N/2
        mlab = mlab - grids.θ.N
    end
    if nlab > grids.ζ.N/2
        nlab = nlab - grids.ζ.N
    end


    #cka does it better than us tbh.
    #max_mode is the index of the maximum mode, essentially useless atm.
    return phi, phi_ft, rgrid[rm[max_mode]], (mlab + grids.θ.pf, nlab + grids.ζ.pf)

end
