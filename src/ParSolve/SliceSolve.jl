"""
    par_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat, solver::SliceSolverT, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec)

Solve the eigenvalue problem Wϕ = ω^2Iϕ using multiple shift and invert transformation to get 'slices' of the spectrum, allowing efficient solving of a larger portion of the spectrum.
"""
function par_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat, solver::SliceSolverT, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec)

    nconv = 0

    evals = ComplexF64[]

    for target in solver.targets
        eps = create_eps(W, I; auto_setup=true)

        EPSSetTarget(eps, target)

        solve!(eps)

        local_conv = EPSGetConverged(eps)

        par_sols_to_file(eps, dir, vecr, veci, local_conv, evals)

        nconv += local_conv
        destroy!(eps)

    end

    #ideally we do this each slice, but that may be causing issues.
    write_evals(evals, dir)
    return nconv

end
