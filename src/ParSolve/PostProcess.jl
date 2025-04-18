"""
    par_sols_to_file(eps::SlepcWrap.SlepcEPS, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32)

Writes the solutions to file. Case for a single shift and invert solve.
"""
function par_sols_to_file(eps::SlepcWrap.SlepcEPS, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        #where the efuncs are stored before they are processed.
        mkpath(dir*"/efuncs_raw") 
    end
    evals = ComplexF64[]

    for ieig in 0:nconv-1

        #the get eigenpair function is kinda weird, seems like it must return stuff
        #hopefully not going to be a huge memory sink
        #in theory we dont need the imaginary part, 
        #should be able to pass null in somehow
        eval, _, vecr, veci = EPSGetEigenpair(eps, ieig, vecr, veci)

        push!(evals, eval)

        efunc_label = "efuncs_raw/"*@sprintf("efunc%05d.hdf5", ieig+1)

        #Creating a viewer each time prevents issues writing to file
        #but is not ideal.
        viewer = PetscViewerHDF5Open(MPI.COMM_WORLD, dir*efunc_label, FILE_MODE_WRITE)

        VecView(vecr, viewer)
        
        PetscViewerDestroy(viewer)

    end

    write_evals(evals, dir)

end

"""
    par_sols_to_file(eps::SlepcWrap.SlepcEPS, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32, evals::Array{ComplexF64}

Writes the solutions to file. Case for slice solve where eigenvalues are added to an array after each solve.
"""
function par_sols_to_file(eps::SlepcWrap.SlepcEPS, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32, evals::Array{ComplexF64})

    for ieig in 0:nconv-1


        #the get eigenpair function is kinda weird, seems like it must return stuff
        #hopefully not going to be a huge memory sink
        #in theory we dont need the imaginary part, 
        #should be able to pass null in somehow
        eval, _, vecr, veci = EPSGetEigenpair(eps, ieig, vecr, veci)

        push!(evals, eval)

        efunc_label = "efuncs_raw/"*@sprintf("efunc%05d.hdf5", length(evals))

        #Creating a viewer each time prevents issues writing to file
        #but is not ideal.
        viewer = PetscViewerHDF5Open(MPI.COMM_WORLD, dir*efunc_label, FILE_MODE_WRITE)

        VecView(vecr, viewer)

        PetscViewerDestroy(viewer)

    end
    #eigenvalues are only written once they have all been found.
    #done in sliceSolve.
end

"""
    par_post_process(dir::String)

Post processes the parallel solutions that have been written to file.
This is done outisde the main functions as there is a bug when converting petsc vec's to julia arrays.
This function is done with a single processor either way.
"""
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

    
    ω = Array{ComplexF64}(undef, nevals)
    mode_labs = Tuple{Int, Int}[] 

    for i in 1:nevals


        efunc_read = @sprintf("efunc%05d.hdf5", i)
        efunc_write = @sprintf("efunc%05d.jld2", i)

        #unfort doesn't handle complex numbers v well
        efunc_split = load_object(dir*"/efuncs_raw/"*efunc_read)

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



"""
    write_evals(evals, dir::String)

Writes the eigenvalues to file.
"""
function write_evals(evals::Array{ComplexF64}, dir::String)

    save_object(dir * "vals_raw.jld2", evals)

end
