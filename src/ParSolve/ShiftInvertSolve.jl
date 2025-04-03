"""
    par_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat)

Solves generalised eigenvalue problem in parallel using Slepc.
"""
function par_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat, solver::ShiftInvertSolverT, dir::String, vecr::PetscWrap.PetscVec, veci::PetscWrap.PetscVec)


    #we create the eps object, auto setup uses the slepcargs above.
    eps = create_eps(W, I; auto_setup=true)

    #this seems to work perfectly!
    #we may still have to create and destroy the eps etc.
    #note target freq should be allowed to be complex, but whatever
    EPSSetTarget(eps, solver.target)

    solve!(eps)



    nconv = EPSGetConverged(eps)

    #vectors are used to store the eigenfunctions
    #vecr, veci = MatCreateVecs(W)

    #perhaps this should be called write to file.
    #no post processing is actually happening.
    par_sols_to_file(eps, dir, vecr, veci, nconv)
    
    destroy!(eps)

    return nconv
end

