
function par_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat, solver::IntervalSolverT, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec)


    eps = create_eps(W, I; auto_setup=true)

    solve!(eps)

    #gets the number of converged solutions
    nconv = EPSGetConverged(eps)


    #write solutions to file, post_processing is done in serial.
    par_sols_to_file(eps, dir, vecr, veci, nconv)
    
    destroy!(eps)

    return nconv
end
