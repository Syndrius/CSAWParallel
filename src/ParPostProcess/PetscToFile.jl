
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
    errs = Float64[]

    tol, maxits = EPSGetTolerances(eps)

    nfound = 0
    for ieig in 0:nconv-1

        #the get eigenpair function is kinda weird, seems like it must return stuff
        #hopefully not going to be a huge memory sink
        #in theory we dont need the imaginary part, 
        #should be able to pass null in somehow
        eval, _, vecr, veci = EPSGetEigenpair(eps, ieig, vecr, veci)

        #think this is actually kind of pointless, as it would require that we 
        #store every eigenfunction as a compairosn.
        #although cka does not compare eigenvlaues, #just seems possible
        #that two branches would have the same eigenvalue? Perhaps they should always bee a wee bit different.
        err = EPSComputeError(eps, ieig)

        #this needs to be done in the postprocessing step if at all.
        #cka compares the equility by taking the dot product of the two and comparing with the norm
        # I think we need to do this is post-processing and just read the eval in,
        #it shouldn't happen that often.
        #ind = new_solution(eps, eval, err, evals, errs, 10 * tol)

        push!(errs, err)
        push!(evals, eval)

        efunc_label = "efuncs_raw/"*@sprintf("efunc%05d.hdf5", ieig+1)

        #Creating a viewer each time prevents issues writing to file
        #but is not ideal.
        viewer = PetscViewerHDF5Open(MPI.COMM_WORLD, dir*efunc_label, FILE_MODE_WRITE)

        VecView(vecr, viewer)
    
        PetscViewerDestroy(viewer)



    end

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        write_evals(evals, dir)
        push!(errs, tol) #just so we have the tolerance used.
        write_errs(errs, dir)
    end

end

"""
    par_sols_to_file(eps::SlepcWrap.SlepcEPS, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32, evals::Array{ComplexF64}

Writes the solutions to file. Case for slice solve where eigenvalues are added to an array after each solve.
"""
function par_sols_to_file(eps::SlepcWrap.SlepcEPS, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec, nconv::Int32, evals::Array{ComplexF64}, errs::Array{Float64})

    #this will get triggered mutiple times.
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        #where the efuncs are stored before they are processed.
        mkpath(dir*"/efuncs_raw") 
    end

    for ieig in 0:nconv-1


        #the get eigenpair function is kinda weird, seems like it must return stuff
        #hopefully not going to be a huge memory sink
        #in theory we dont need the imaginary part, 
        #should be able to pass null in somehow
        eval, _, vecr, veci = EPSGetEigenpair(eps, ieig, vecr, veci)

        err = EPSComputeError(eps, ieig)

        push!(errs, err)

        push!(evals, eval)

        efunc_label = "efuncs_raw/"*@sprintf("efunc%05d.hdf5", length(evals))

        #Creating a viewer each time prevents issues writing to file
        #but is not ideal.
        viewer = PetscViewerHDF5Open(MPI.COMM_WORLD, dir*efunc_label, FILE_MODE_WRITE)

        VecView(vecr, viewer)

        PetscViewerDestroy(viewer)

    end
    #eigenvalues are written multiple times, so that if errors occur, some of the slices are still saved.
    #write_evals(evals, dir)
end





"""
    write_evals(evals, dir::String)

Writes the eigenvalues to file.
"""
function write_evals(evals::Array{ComplexF64}, dir::String)

    save_object(dir * "vals_raw.jld2", evals)

end

function write_errs(errs::Array{Float64}, dir::String)

    save_object(dir * "errs.jld2", errs)

end


#this might need to be moved!
function EPSComputeError(eps::SlepcWrap.SlepcEPS, ind::Int64)
    rel_err = Ref{PetscReal}(0)


    #1 here is the EPS_ERROR_RELATIVE
    error = ccall((:EPSComputeError, SlepcWrap.libslepc), PetscErrorCode, (Ptr{Cvoid}, PetscInt, PetscInt, Ref{PetscReal}), eps, PetscInt(ind), PetscInt(1), rel_err)
    @assert iszero(error)

    return rel_err[]
end


