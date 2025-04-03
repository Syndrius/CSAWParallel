
#this is now working at a basic level, could still do with extra stuff
#eg check for double ups and take lowest errer efunc etc. For now we can probably just assume the targets are chosen well

#function that slicing up a region and solves with shift and invert multiple times
#to build up a larger portion of the spectrum

function par_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat, solver::SliceSolverT, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec)

    nconv = 0

    evals = ComplexF64[]

    #think we will just create and destroy the eps object everytime, not efficient, but defs will work.
    for target in solver.targets
        eps = create_eps(W, I; auto_setup=true)

        EPSSetTarget(eps, target)

        solve!(eps)

        local_conv = EPSGetConverged(eps)


        par_sols_to_file(eps, dir, vecr, veci, local_conv, evals)

        nconv += local_conv
        destroy!(eps)

    end

    write_evals(evals, dir)
    return nconv

end
