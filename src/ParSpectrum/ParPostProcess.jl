
"""
    function par_post_process(eps::SlepcWrap.SlepcEPS, dir::String, grids::MID.GridsT, geo::MID.GeoParamsT, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32)

Processes the eigenvalues and eigenfunctions. This includes normalising the eigenvalues, converting the eigenfunctions back to the 3d grid, computing the continuum reconstruction and writing the outputs to file.
"""
function par_post_process(eps::SlepcWrap.SlepcEPS, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        #where the efuncs are stored before they are processed.
        mkpath(dir*"/efuncs_raw") 
    end

    #viewer = PetscViewerCreate()
    


    for ieig in 0:nconv-1

        #the extra things returned here are the complex parts,
        #but in scalar_type=complex, the real parts are also comlpex
        #and complex parts are just zero.
        #this is returning the parallel eigenvalue and the parallel eigenvector
        #in Petsc format.
        #changing the inputs for this made no difference.
        EPSGetEigenvector(eps, ieig, vecr, veci)

        efunc_label = "efuncs_raw/"*@sprintf("efunc%04d.hdf5", ieig+1)

        #certainly not ideal to be creating a viewer each time..
        viewer = PetscViewerHDF5Open(MPI.COMM_WORLD, dir*efunc_label, FILE_MODE_WRITE)
        #not certain if we want to write to different files
        #may want to just write to a single file but different group
        #could be worth trying both
        #PetscViewerFileSetMode(viewer) #defaults to FILE_MODE_WRITE
        #PetscViewerSetType(viewer, "hdf5")
        #PetscViewerFileSetName(viewer, dir*efunc_label)

        VecView(vecr, viewer)

    end
    #PetscViewerDestroy(viewer)

end

#same function but the option of a deriv.
#probably will replace the old version soon.
function process_hdf5(dir::String, deriv::Bool=false)

    #if isfile(dir*"/evals.jld2")
        #this will not be good enough when we also want to find a derivative!
        #display("Processing already complete.")
        #return
    #end

    #these almost certianly already exist, but in case they don't
    if deriv
        mkpath(dir*"/efuncs_deriv")
        mkpath(dir*"/efuncs_ft_deriv")
    else
        mkpath(dir*"/efuncs")
        mkpath(dir*"/efuncs_ft")
    end
        
    #mkpath(dir*"/efuncs_raw")

    #this will work for the newest cases but not the old ones.
    #probably need a try catch??
    prob, grids = inputs_from_file(dir=dir)


    #maybe better to do this via hdf5 or something.
    vals = open(dir*"vals.dat", "r") do file
        s = readlines(file)[2:end-1] 
        parse.(ComplexF64, s)
    end

    #number of evals/efuncs.
    #nevals = length(vals)

    #this accounts for cases when no every single efuncs get written for some reason.
    strs = readdir(dir*"/efuncs_raw")
    #not sure if a lock hdf5 is always written?
    #sometimes loc sometimes lock...
    #may be best to just delete the appropriate number of evals
    #or delete the loc vars from efuncs_raw if this happens again.
    #these cases may just be cooked though lol.
    if strs[end][end-2:end] == "loc"

        nevals = Int64(length(strs)/2)
    else
        nevals = length(strs)
    end

    display(nevals)
    #allocates placeholder arrays used during computation.
    ϕp, ϕpft = MID.PostProcessing.allocate_phi_arrays(grids, deriv=deriv)

    rms = Array{Float64}(undef, nevals)

    plan = MID.PostProcessing.create_ft_plan(ϕpft, grids)

    rgrid = inst_grids(grids)[1]

    

    #arrays to store the maximum value of ϕft and the correspondning r value.
    rmarray = Array{Int64}(undef, grids.θ.N, grids.ζ.N)
    ϕmarray = Array{Float64}(undef, grids.θ.N, grids.ζ.N)

    
    #rmode = zeros(Float64, nconv)
    #mlab = zeros(Int64, nconv)
    #nlab = zeros(Int64, nconv)
    ω = Array{ComplexF64}(undef, nevals)
    mode_labs = Tuple{Int, Int}[] 

    for i in 1:nevals


        #this is fkn slow af.
        #need to stick with jld2 I think.
        efunc_read = @sprintf("efunc%04d.hdf5", i)
        efunc_write = @sprintf("efunc%04d.jld2", i)

        #unfort doesn't handle complex numbers v well
        efunc_split = load_object(dir*"/efuncs_raw/"*efunc_read)#[1, :]

        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im


        MID.PostProcessing.reconstruct_phi!(efunc, grids, ϕp, ϕpft, plan)
        
        rind, mode_lab = MID.PostProcessing.label_mode(ϕpft, grids, rmarray, ϕmarray)


        #if deriv
        #    rm, modelab = MID.Spectrum.label_mode(ϕft[:, :, :, 1], grids)
        #else
        #    rm, modelab = MID.Spectrum.label_mode(ϕft, grids)
        #end

        
        #phi, phi_ft, rm, modelab = process_efunc(efunc, grids)

        push!(mode_labs, mode_lab)

        rms[i] = rgrid[rind]

        #normalise the eigenvalues
        ω[i] = prob.geo.R0 * sqrt(vals[i])
        if deriv
            save_object(dir * "/efuncs_deriv/"*efunc_write, ϕp)
            save_object(dir * "/efuncs_ft_deriv/"*efunc_write, ϕpft)
        else
            save_object(dir * "/efuncs/"*efunc_write, ϕp)
            save_object(dir * "/efuncs_ft/"*efunc_write, ϕpft)
        end
        
    end
    evals = MID.EvalsT(ω, rms, mode_labs)

    save_object(dir*"/evals.jld2", evals)

end



#used to process the direct output of slepc for when the weird stuff happens.
#this should be called in serial!
#maybe we should ignore jld2 and just use hdf5? but jld2 is nicer...
#so this didn't seem to work...
#unclear why tbh.
#don't understand what has changed, but this is taking way longer than before....
function process_hdf5_old(dir::String)

    if isfile(dir*"/evals.jld2")
        display("Processing already complete.")
        return
    end

    #these almost certianly already exist, but in case they don't
    mkpath(dir*"/efuncs")
    mkpath(dir*"/efuncs_ft")
    #mkpath(dir*"/efuncs_raw")

    #this will work for the newest cases but not the old ones.
    #probably need a try catch??
    prob, grids = inputs_from_file(dir=dir)


    #maybe better to do this via hdf5 or something.
    vals = open(dir*"vals.dat", "r") do file
        s = readlines(file)[2:end-1] 
        parse.(ComplexF64, s)
    end


    nconv = length(vals)

    #ideally we won't do this.
    #evals_old = load(dir*"evals.jld2")
    #nconv = length(evals.ω)

    


    #ie eval, rmax, label??
    #label isn't a float though...
    #cont_reconstruction = zeros(Float64, nconv, nconv, nconv)
    rmode = zeros(Float64, nconv)
    #mlab = zeros(Int64, nconv)
    #nlab = zeros(Int64, nconv)
    ω = zeros(ComplexF64, nconv)
    modelabs = Tuple{Int, Int}[] 

    for i in 1:nconv

        #this is fkn slow af.
        #need to stick with jld2 I think.
        efunc_read = @sprintf("efunc%04d.hdf5", i)
        efunc_write = @sprintf("efunc%04d.jld2", i)

        #this gives a warning but seems to be fine
        #this is much faster, but still takes sometime, 
        #may want to do this as a job and parralise?
        #could be a bu=it extra though lol
        #not sure why this was so much faster before?
        #perhaps we had a better node?
        #we can probably optimise this a bit by preallocating data and creating fourier transform plans etc.
        #perhaps 2 is the ind we want?
        #why the fk is there two????? golly gosh this is fkn stoopid.
        efunc_split = load_object(dir*"/efuncs_raw/"*efunc_read)#[1, :]

        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im
        #fid = h5open(dir*"/efuncs_raw/"*efunc_read, "r")
        #key = (keys(fid))[1]
        #efunc = read(fid[key])[1, :]
        #close(fid)
        #display(efunc)
        #efunc = load_object(dir*"/efuncs_raw/"*efunc_label)
        phi, phi_ft, rm, modelab = process_efunc(efunc, grids)

        push!(modelabs, modelab)

        rmode[i] = rm

        #normalise the eigenvalues
        ω[i] = prob.geo.R0 * sqrt(vals[i])

        save_object(dir * "/efuncs/"*efunc_write, phi)
        save_object(dir * "/efuncs_ft/"*efunc_write, phi_ft)
        
    end
    evals = MID.EvalsT(ω, rmode, modelabs)

    save_object(dir*"/evals.jld2", evals)

end





function process_efunc(efunc, grids::MID.FSSGridsT)


    
    phi_ft = zeros(ComplexF64, grids.r.N, grids.θ.count, grids.ζ.count)

    θgrid_size = compute_ifft_grid(grids.θ)
    ζgrid_size = compute_ifft_grid(grids.ζ)

    phi = zeros(ComplexF64, grids.r.N, θgrid_size, ζgrid_size)


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

    θgrid = LinRange(0, 2π, θgrid_size+1)[1:end-1]
    ζgrid = LinRange(0, 2π, ζgrid_size+1)[1:end-1]

    #ift but for only the modes specified.
    for j in 1:grids.r.N
        #for k in 1:grids.ζ.count
        for k in 1:1:length(nlist)
            for l in 1:1:length(mlist)

                phi[j, :, :] += phi_ft[j, l, k] .* exp.(1im * mlist[l] .* θgrid .+ 1im * nlist[k] .* ζgrid' )
            end
        end
    end


    #for i in 1:grids.r.N
        #ifft to restore original ϕ
    #    phi[i, :, :] = ifft(phi_ft[i, :, :], [1, 2])
    #end 
    #phi[:, :, :] = ifft(phi_ft[:, :, :], [2, 3])

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

    #phi_ffs = MID.reconstruct_phi(efunc, grids)

    phi_ffs = zeros(ComplexF64, grids.r.N, grids.θ.N, grids.ζ.count)
    phi_ft = zeros(ComplexF64, grids.r.N, grids.θ.N, grids.ζ.count)

    ζgrid_size = compute_ifft_grid(grids.ζ)
    phi = zeros(ComplexF64, grids.r.N, grids.θ.N, ζgrid_size)


    rgrid, θgrid, _, nlist, _= instantiate_grids(grids)

    #maxphi = -100
    #we will probably want a plan for fourier transform.

    m = grids.θ.pf

    #think this can be done more efficeintly by knowing the pattern we should be able to set multiple values at once.
    #but this can do for now.
    #eg have a look at cka's way of doing it.
    for i in 1:4:matrix_dim(grids)

        #note these are the indicies.
        r, θ, ζ, hs = index_to_grid(i, grids)
        phi_ffs[r, θ, ζ] = efunc[i] * exp(1im * m * θgrid[θ])

        #shouldn't need this if we skip over each 8.
        #if hs == 1
            #may be the wrong way around!
            #this doesn't seem to have worked as expected tbh!
             #.* exp(1im * (m * θgrid[θ] + n * ζgrid[ζ]))
        #end
    end




    phi_ft = fft(phi_ffs, [2])
    
    ζgrid = LinRange(0, 2π, ζgrid_size+1)[1:end-1]
    

    for j in 1:grids.r.N
            
        for k in 1:grids.θ.N
            for l in 1:1:length(nlist)

                phi[j, k, :] += phi_ffs[j, k, l] .* exp.(1im * nlist[l] .* ζgrid)
            end
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


    rgrid, θgrid, ζgrid = MID.instantiate_grids(grids)

    #maxphi = -100
    #we will probably want a plan for fourier transform.

    m = grids.θ.pf
    n = grids.ζ.pf

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
            phi[r, θ, ζ] = efunc[i] * exp(1im * (m * θgrid[θ] + n * ζgrid[ζ]))
        end
    end


    phi_ft = fft(phi, [2, 3])


    rm = zeros(Int64, grids.θ.N, grids.ζ.N)
    ϕm = zeros(Float64, grids.θ.N, grids.ζ.N)

    for j in 1:grids.θ.N, k in 1:grids.ζ.N

        rm[j, k] = argmax(abs.(real.(phi_ft[:, j, k])))
        ϕm[j, k] = abs.(real.(phi_ft[rm[j, k], j, k]))
    end

    max_mode = argmax(ϕm) #still need to convert this to the correct label... 


    mlab = MID.mode_label(max_mode[1], grids.θ)
    nlab = MID.mode_label(max_mode[2], grids.ζ)


    #cka does it better than us tbh.
    #max_mode is the index of the maximum mode, essentially useless atm.
    return phi, phi_ft, rgrid[rm[max_mode]], (mlab, nlab)

end
#





"""
    function par_post_process(eps::SlepcWrap.SlepcEPS, dir::String, grids::MID.GridsT, geo::MID.GeoParamsT, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32)

Processes the eigenvalues and eigenfunctions. This includes normalising the eigenvalues, converting the eigenfunctions back to the 3d grid, computing the continuum reconstruction and writing the outputs to file.
"""
function old_par_post_process(eps::SlepcWrap.SlepcEPS, dir::String, grids::MID.GridsT, geo::MID.GeoParamsT, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32)
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
        #mkpath(dir*"/efuncs_raw")


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
        #changing the inputs for this made no difference.
        vp, _, vecp, _ = EPSGetEigenpair(eps, ieig, vecr, veci)

        #display(vpi)
        #this still writes in a stupid formate.
        #viewer = PetscViewerASCIIOpen(MPI.COMM_WORLD, "test.dat")
        #VecView(vecpr, viewer)

        #so this gets a julia array for each proc...
        #create a julia vector for each procs portion of the total vector
        #something wrong with this :(
        #gridap version of this is slightly different. 
        #maybe the solution, will take some effort to figure out though...
        #it is odd that it is so consistent tbh.
        #I guess Gadi's queue structure probably uses the same memory blocks etc
        #maybe with true as the second arg this will work,
        #but it takes at least 4 times longer. so not practical...
        #we will just have to try and use GridAp version.
        #jvecp, ref = VecGetArray(vecp, false)
        #MPI.Barrier(MPI.COMM_WORLD)
        #vecr or vecp here doesn't matter!
        jvecp = convert_to_jl(vecp)
        #MPI.Barrier(MPI.COMM_WORLD)


        #this solution will work, pretty fkn annoying though
        #we will have to postprocess later.
        #viewer = PetscViewerHDF5Open(MPI.COMM_WORLD, dir*"test.hdf5", FILE_MODE_WRITE)

        #VecView(vecp, viewer)




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

            #maybe also try saving normally, i.e. not using jld2
            #save_object(dir * "/efuncs_raw/"*efunc_label, jvec)

            #cont_reconstruction[i]

            #while not ideal, I think this is the way forward, now fourier transform etc and write to file
            #ideally want to get phi structure in file form. 
            #display(jvec)
            #save_object("test.jld2", g_jvec)
        end

        
        #VecRestoreArray(vecp, ref)
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

#hopefully fix the vecgetarray problem.
#this did not work.
function convert_to_jl(vec::PetscVec)

    r_pv = Ref{Ptr{PetscScalar}}()
    error = ccall((:VecGetArray, PetscWrap.libpetsc), PetscErrorCode, (CVec, Ref{Ptr{PetscScalar}}), vec, r_pv)

    @assert iszero(error)

    r_sz = Ref{PetscInt}()
    
    error = ccall((:VecGetLocalSize, PetscWrap.libpetsc), PetscErrorCode, (CVec, Ref{PetscInt}), vec, r_sz)

    @assert iszero(error)

    v = unsafe_wrap(Array, r_pv[], r_sz[]; own=false)

    error = ccall((:VecRestoreArray, PetscWrap.libpetsc), PetscErrorCode, (CVec, Ptr{Ptr{PetscScalar}}), vec, Ref(pointer(v)))

    @assert iszero(error)
    return v
end
