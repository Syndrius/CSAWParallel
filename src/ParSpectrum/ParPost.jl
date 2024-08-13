
#file for postprocessing, eg fourier transform and continuum etc.

#may just combine this into solve, otherwise solve does fk all.

#no fkn idea what the type is lol
function ParPost(eps)
    #unclear if this can be done in parallel or not.
    #looks like it may be possible???
    #think start by just getting root to do it all, then maybe we want to split up the solutions.
    

    #will need different versions of this for different grids.
    #looks like a single core is the way.
    #may be possible to parallise later but may not be worth it.

    #number of solutions should be comparable to work of single core anyway.

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        nconv = EPSGetConverged(eps)
        display("Number of converged eigenvalues = ", nconv)
    end

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0

        for i in 0:nconv-1

        
    end







end