#this is the case for a single target
function par_sols_to_file(eps::SlepcWrap.SlepcEPS, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        #where the efuncs are stored before they are processed.
        mkpath(dir*"/efuncs_raw") 
    end
    evals = ComplexF64[]

    display("This file is being called")

    for ieig in 0:nconv-1

        #the extra things returned here are the complex parts,
        #but in scalar_type=complex, the real parts are also comlpex
        #and complex parts are just zero.
        #this is returning the parallel eigenvalue and the parallel eigenvector
        #in Petsc format.
        #changing the inputs for this made no difference.
        #EPSGetEigenvector(eps, ieig, vecr, veci)

        #the get eigenpair function is kinda weird, seems like it must return stuff
        #hopefully not going to be a huge memory sink
        #in theory we dont need the imaginary part, 
        #should be able to pass null in somehow
        eval, _, vecr, veci = EPSGetEigenpair(eps, ieig, vecr, veci)

        push!(evals, eval)

        efunc_label = "efuncs_raw/"*@sprintf("efunc%05d.hdf5", ieig+1)

        #certainly not ideal to be creating a viewer each time..
        viewer = PetscViewerHDF5Open(MPI.COMM_WORLD, dir*efunc_label, FILE_MODE_WRITE)
        #viewer = PetscViewerASCIIOpen(MPI.COMM_WORLD, dir*efunc_label)
        #not certain if we want to write to different files
        #may want to just write to a single file but different group
        #could be worth trying both
        #PetscViewerFileSetMode(viewer) #defaults to FILE_MODE_WRITE
        #PetscViewerSetType(viewer, "hdf5")
        #PetscViewerFileSetName(viewer, dir*efunc_label)

        VecView(vecr, viewer)
        #maybe the viewer is running our of memory???
        #given we are able to solve but not write them all??
        #and all the previous ones still have the .loc thing?
        #this may have made a big difference for memory usage...
        PetscViewerDestroy(viewer)

    end

    write_evals(evals, dir)

end

#multiple dispatch for the different solver methods is simplest I think.
function par_sols_to_file(eps::SlepcWrap.SlepcEPS, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32, evals::Array{ComplexF64})

    for ieig in 0:nconv-1

        #display(ieig)

        #the extra things returned here are the complex parts,
        #but in scalar_type=complex, the real parts are also comlpex
        #and complex parts are just zero.
        #this is returning the parallel eigenvalue and the parallel eigenvector
        #in Petsc format.
        #changing the inputs for this made no difference.
        #EPSGetEigenvector(eps, ieig, vecr, veci)

        #the get eigenpair function is kinda weird, seems like it must return stuff
        #hopefully not going to be a huge memory sink
        eval, _, vecr, veci = EPSGetEigenpair(eps, ieig, vecr, veci)

        push!(evals, eval)

        efunc_label = "efuncs_raw/"*@sprintf("efunc%05d.hdf5", length(evals))

        #certainly not ideal to be creating a viewer each time..
        viewer = PetscViewerHDF5Open(MPI.COMM_WORLD, dir*efunc_label, FILE_MODE_WRITE)
        #viewer = PetscViewerASCIIOpen(MPI.COMM_WORLD, dir*efunc_label)
        #not certain if we want to write to different files
        #may want to just write to a single file but different group
        #could be worth trying both
        #PetscViewerFileSetMode(viewer) #defaults to FILE_MODE_WRITE
        #PetscViewerSetType(viewer, "hdf5")
        #PetscViewerFileSetName(viewer, dir*efunc_label)

        VecView(vecr, viewer)
        #maybe the viewer is running our of memory???
        #given we are able to solve but not write them all??
        #and all the previous ones still have the .loc thing?
        #this may have made a big difference for memory usage...
        PetscViewerDestroy(viewer)

    end
    #eigenvalues are only written once they have all been found.
end

#this is process hdf5
#ignoring the derivative stuff for now!
#note that this function needs to be called in serial
function par_post_process(dir::String)
    mkpath(dir*"/efuncs")
    mkpath(dir*"/efuncs_ft")

    prob, grids, _ = inputs_from_file(dir=dir)

    vals = load_object(dir*"vals_raw.jld2")
    nevals = length(vals)

    ϕp, ϕpft = PostProcessing.allocate_phi_arrays(grids, deriv=false)
    rms = Array{Float64}(undef, nevals)

    plan = PostProcessing.create_ft_plan(ϕpft, grids)

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
        efunc_read = @sprintf("efunc%05d.hdf5", i)
        efunc_write = @sprintf("efunc%05d.jld2", i)

        #unfort doesn't handle complex numbers v well
        efunc_split = load_object(dir*"/efuncs_raw/"*efunc_read)#[1, :]

        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im


        PostProcessing.reconstruct_phi!(efunc, grids, ϕp, ϕpft, plan)
        
        rind, mode_lab = PostProcessing.label_mode(ϕpft, grids, rmarray, ϕmarray)


        push!(mode_labs, mode_lab)

        rms[i] = rgrid[rind]

        #normalise the eigenvalues
        ω[i] = prob.geo.R0 * sqrt(vals[i])

        save_object(dir * "/efuncs/"*efunc_write, ϕp)
        save_object(dir * "/efuncs_ft/"*efunc_write, ϕpft)
        
    end
    evals = EvalsT(ω, rms, mode_labs)

    save_object(dir*"/evals.jld2", evals)
end



function write_evals(evals, dir::String)

    save_object(dir * "vals_raw.jld2", evals)

end

#######################################################################
#older versions

"""
    function par_post_process(eps::SlepcWrap.SlepcEPS, dir::String, grids::MID.GridsT, geo::MID.GeoParamsT, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32)

Processes the eigenvalues and eigenfunctions. This includes normalising the eigenvalues, converting the eigenfunctions back to the 3d grid, computing the continuum reconstruction and writing the outputs to file.
"""
function old_par_post_process(eps::SlepcWrap.SlepcEPS, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32, evals::Array{ComplexF64}=ComplexF64[])

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        #where the efuncs are stored before they are processed.
        mkpath(dir*"/efuncs_raw") 
    end

    #viewer = PetscViewerCreate()
    slice = false
    
    #testing, the condition needs ot be better
    #if !isnothing(evals)
    if !slice

        evals = ComplexF64[]

        for ieig in 0:nconv-1

            #the extra things returned here are the complex parts,
            #but in scalar_type=complex, the real parts are also comlpex
            #and complex parts are just zero.
            #this is returning the parallel eigenvalue and the parallel eigenvector
            #in Petsc format.
            #changing the inputs for this made no difference.
            #EPSGetEigenvector(eps, ieig, vecr, veci)

            #the get eigenpair function is kinda weird, seems like it must return stuff
            #hopefully not going to be a huge memory sink
            #in theory we dont need the imaginary part, 
            #should be able to pass null in somehow
            eval, _, vecr, veci = EPSGetEigenpair(eps, ieig, vecr, veci)

            push!(evals, eval)

            efunc_label = "efuncs_raw/"*@sprintf("efunc%05d.hdf5", ieig+1)

            #certainly not ideal to be creating a viewer each time..
            viewer = PetscViewerHDF5Open(MPI.COMM_WORLD, dir*efunc_label, FILE_MODE_WRITE)
            #viewer = PetscViewerASCIIOpen(MPI.COMM_WORLD, dir*efunc_label)
            #not certain if we want to write to different files
            #may want to just write to a single file but different group
            #could be worth trying both
            #PetscViewerFileSetMode(viewer) #defaults to FILE_MODE_WRITE
            #PetscViewerSetType(viewer, "hdf5")
            #PetscViewerFileSetName(viewer, dir*efunc_label)

            VecView(vecr, viewer)
            #maybe the viewer is running our of memory???
            #given we are able to solve but not write them all??
            #and all the previous ones still have the .loc thing?
            #this may have made a big difference for memory usage...
            PetscViewerDestroy(viewer)

        end

        write_evals(evals, dir)
    else

        for ieig in 0:nconv-1

            #display(ieig)

            #the extra things returned here are the complex parts,
            #but in scalar_type=complex, the real parts are also comlpex
            #and complex parts are just zero.
            #this is returning the parallel eigenvalue and the parallel eigenvector
            #in Petsc format.
            #changing the inputs for this made no difference.
            #EPSGetEigenvector(eps, ieig, vecr, veci)

            #the get eigenpair function is kinda weird, seems like it must return stuff
            #hopefully not going to be a huge memory sink
            eval, _, vecr, veci = EPSGetEigenpair(eps, ieig, vecr, veci)

            push!(evals, eval)

            efunc_label = "efuncs_raw/"*@sprintf("efunc%05d.hdf5", length(evals))

            #certainly not ideal to be creating a viewer each time..
            viewer = PetscViewerHDF5Open(MPI.COMM_WORLD, dir*efunc_label, FILE_MODE_WRITE)
            #viewer = PetscViewerASCIIOpen(MPI.COMM_WORLD, dir*efunc_label)
            #not certain if we want to write to different files
            #may want to just write to a single file but different group
            #could be worth trying both
            #PetscViewerFileSetMode(viewer) #defaults to FILE_MODE_WRITE
            #PetscViewerSetType(viewer, "hdf5")
            #PetscViewerFileSetName(viewer, dir*efunc_label)

            VecView(vecr, viewer)
            #maybe the viewer is running our of memory???
            #given we are able to solve but not write them all??
            #and all the previous ones still have the .loc thing?
            #this may have made a big difference for memory usage...
            PetscViewerDestroy(viewer)

        end
        
    end
    #PetscViewerDestroy(viewer)

end


"""
    function par_post_process(eps::SlepcWrap.SlepcEPS, dir::String, grids::MID.GridsT, geo::MID.GeoParamsT, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32)

Processes the eigenvalues and eigenfunctions. This includes normalising the eigenvalues, converting the eigenfunctions back to the 3d grid, computing the continuum reconstruction and writing the outputs to file.
"""
function oldest_par_post_process(eps::SlepcWrap.SlepcEPS, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32)

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

        efunc_label = "efuncs_raw/"*@sprintf("efunc%05d.hdf5", ieig+1)

        #certainly not ideal to be creating a viewer each time..
        viewer = PetscViewerHDF5Open(MPI.COMM_WORLD, dir*efunc_label, FILE_MODE_WRITE)
        #viewer = PetscViewerASCIIOpen(MPI.COMM_WORLD, dir*efunc_label)
        #not certain if we want to write to different files
        #may want to just write to a single file but different group
        #could be worth trying both
        #PetscViewerFileSetMode(viewer) #defaults to FILE_MODE_WRITE
        #PetscViewerSetType(viewer, "hdf5")
        #PetscViewerFileSetName(viewer, dir*efunc_label)

        VecView(vecr, viewer)
        #maybe the viewer is running our of memory???
        #given we are able to solve but not write them all??
        #and all the previous ones still have the .loc thing?
        #this may have made a big difference for memory usage...
        PetscViewerDestroy(viewer)

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
    #vals = open(dir*"vals.dat", "r") do file
    #    s = readlines(file)[2:end-1] 
    #    parse.(ComplexF64, s)
    #end
    vals = load_object(dir*"vals_raw.jld2")

    #number of evals/efuncs.
    #nevals = length(vals)

    #this accounts for cases when no every single efuncs get written for some reason.
    strs = readdir(dir*"efuncs_raw")
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
    nevals = length(vals)
    #allocates placeholder arrays used during computation.
    #not idea to import like this tbh!, but Postprocessing doesn't export any of this!
    ϕp, ϕpft = PostProcessing.allocate_phi_arrays(grids, deriv=deriv)

    rms = Array{Float64}(undef, nevals)

    plan = PostProcessing.create_ft_plan(ϕpft, grids)

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
        efunc_read = @sprintf("efunc%05d.hdf5", i)
        efunc_write = @sprintf("efunc%05d.jld2", i)

        #unfort doesn't handle complex numbers v well
        efunc_split = load_object(dir*"/efuncs_raw/"*efunc_read)#[1, :]

        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im


        PostProcessing.reconstruct_phi!(efunc, grids, ϕp, ϕpft, plan)
        
        rind, mode_lab = PostProcessing.label_mode(ϕpft, grids, rmarray, ϕmarray)


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
    evals = EvalsT(ω, rms, mode_labs)

    save_object(dir*"/evals.jld2", evals)

end



