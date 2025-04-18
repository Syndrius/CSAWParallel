"""
    par_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat, solver::ShiftInvertSolverT, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec)

Solves generalised eigenvalue problem in parallel using Slepc by using a shift and invert transformation to find the eigenvalus nearest to the target.
"""
function par_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat, solver::ShiftInvertSolverT, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec)


    #we create the eps object, auto setup uses the slepcargs 
    eps = create_eps(W, I; auto_setup=true)

    EPSSetTarget(eps, solver.target)

    solve!(eps)

    #gets the number of converged solutions
    nconv = EPSGetConverged(eps)


    #write solutions to file, post_processing is done in serial.
    par_sols_to_file(eps, dir, vecr, veci, nconv)
    
    destroy!(eps)

    return nconv
end

