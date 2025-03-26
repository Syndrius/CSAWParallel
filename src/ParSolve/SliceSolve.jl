
#TODO

#function that slicing up a region and solves with shift and invert multiple times
#to build up a larger portion of the spectrum
function slice_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat, targets::Array{Float64}, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec)

    nconv = 0

    evals = ComplexF64[]

    #think we will just create and destroy the eps object everytime, not efficient, but defs will work.
    for target in targets
        eps = create_eps(W, I; auto_setup=true)

        EPSSetTarget(eps, target)

        solve!(eps)

        local_conv = EPSGetConverged(eps)


        par_post_process(eps, dir, vecr, veci, local_conv, evals)

        nconv += local_conv
        destroy!(eps)

    end

    write_evals(evals, dir)
    return nconv

end
