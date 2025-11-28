
"""
    solve(P::PetscWrap.PetscMat, Q::PetscWrap.PetscMat, solver::IntervalSolverT, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec)

Solves the eigenvalue problem by finding all solutions in between the given interval.
Work in progress...
Slice solve is more effective.
"""
function solve(P::PetscWrap.PetscMat, Q::PetscWrap.PetscMat, solver::IntervalSolverT, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec)

    eps = create_eps(P, Q; auto_setup=true)

    solve!(eps)

    #gets the number of converged solutions
    nconv = EPSGetConverged(eps)

    #write solutions to file, post_processing is done in serial.
    par_sols_to_file(eps, dir, vecr, veci, nconv)
    
    destroy!(eps)

    return nconv
end
