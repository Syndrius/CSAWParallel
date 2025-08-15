"""
    par_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat, solver::SliceSolverT, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec)

Solve the eigenvalue problem Wϕ = ω^2Iϕ using multiple shift and invert transformation to get 'slices' of the spectrum, allowing efficient solving of a larger portion of the spectrum.
"""
function par_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat, solver::SliceSolverT, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec)

    nconv = 0

    evals = ComplexF64[]
    errs = Float64[]

    tol = 0.0
    for target in solver.targets
        eps = create_eps(W, I; auto_setup=true)

        EPSSetTarget(eps, target)

        solve!(eps)

        local_conv = EPSGetConverged(eps)

        par_sols_to_file(eps, dir, vecr, veci, local_conv, evals, errs)

        nconv += local_conv
        tol, maxits = EPSGetTolerances(eps)
        destroy!(eps)
    end

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        write_evals(evals, dir)
        push!(errs, tol) #just so we have the tolerance used.
        write_errs(errs, dir)
    end

    return nconv

end
