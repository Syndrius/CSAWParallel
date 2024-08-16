"""
    par_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat)

Solves generalised eigenvalue problem in parallel using Slepc.
"""
function par_solve(W::PetscWrap.PetscMat, I::PetscWrap.PetscMat)

    #we create the eps object, auto setup uses the slepcargs above.
    eps = create_eps(W, I; auto_setup=true)

    solve!(eps)

    return eps
end

