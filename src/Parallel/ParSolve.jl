
#solves in parallel using Slepc

function par_solve(rows, cols, Wdata, Idata; σ, nev=5, R0, n, dir::String)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    root = 0
    
    #display("Big Boi")
    #eps_nev doesn't seem to do anything???
    slepcargs = @sprintf("-eps_target %s -eps_nev %s -st_pc_factor_shift_type NONZERO -st_type sinvert ", σ, nev)
    #running in parallel with -st_type sinvert gives very large evals, with eps_harmonic we get ve small evals, neither is good!
    #slepcargs = @sprintf("-eps_nev 3 -eps_target %s -eps_view -st_pc_factor_shift_type NONZERO", σ)
    
    #if rank == root
    #    display(slepcargs)
    #end
    SlepcInitialize(slepcargs)
    #SlepcInitialize("-eps_nev 3")# -eps_type lobpcg")
    #n = 2*N * grids.pmd.count * grids.tmd.count
    #n = 4 * 1 * 2

    W = create_matrix(n, n, auto_setup=true)
    I = create_matrix(n, n, auto_setup=true)

    
    set_values!(W, cols, rows, Wdata)
    set_values!(I, cols, rows, Idata)
    
    #MatSetValues(W, length_data, slep_rows, length_data, slep_cols, slep_Wdata, INSERT_VALUES);

    assemble!(W)
    assemble!(I)

    #display(MatGetLocalSize(W))
    display(MatGetOwnershipRangeColumn(W))
    
    #MatView(I)
    #
    #MPI.Barrier(comm)
    eps = create_eps(W, I; auto_setup=true)

    #so this doesn't seem to work on Pc now...
    solve!(eps)
    

    nconv = neigs(eps)

    # Then we can get/display these eigenvalues (more precisely their square root, i.e ``\simeq \omega``)
    #note nconv != nev, ie more can converge!
    #if rank==root
    for ieig in 1:nconv
        eig = get_eig(eps, ieig)
        display(eig)
    end
    #end

    # Export eigenvalues to a file
    #not sure that the eigenvectors is working as indended, but evals are working well!
    eigenvalues2file(eps, dir*"eigvals.dat", 1:nconv)

    # Export eigenvectors as two ASCII matrices (real/imag) (experimental function)
    eigenvectors2file(eps, dir*"eigfuncs.dat", 1:nconv)
    
    #display(10.0 .* sqrt.(get_eigenvalues(eps)))


    destroy!(W)
    destroy!(I)
    destroy!(eps)

    SlepcFinalize()
    


end

