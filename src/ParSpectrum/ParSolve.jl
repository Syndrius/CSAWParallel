"""
Solves generalised eigenvalue problem in parallel using Slepc.

# Args
- rows::Vector{Int64} Array containing the row indexes.
- cols::Vector{Int64} Array containing the column indexes.
- Wdata::Vector{ComplexF64} Array containing the data for W, corresponding to the rows and cols.
- Idata::Vector{ComplexF64} Array containing the data for I, corresponding to the rows and cols.
- σ::Float64 Frequency for 'shift and invert', typically this is the TAE frequency.
- nev::Int64 Number of eigenvalues to solve for.
- n::Int64 Size of the matrix.
- dir::String Directory the results are written to.
"""
function par_solve(W, I)


    
    
    #matrices are then assembled.
    assemble!(W)
    assemble!(I)

    #we create the eps object, auto setup uses the slepcargs above.
    eps = create_eps(W, I; auto_setup=true)

    
    #solve the problem.
    #slepc args above will automatically write the solution to file.
    #note sometimes nev is ignored and solution just gives us the number of converged eigenvalues.
    solve!(eps)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display(EPSGetConverged(eps))
    end
    
    #free the memory used.
    
    destroy!(eps)

    
end


"""
Solves generalised eigenvalue problem in parallel using Slepc.

# Args
- rows::Vector{Int64} Array containing the row indexes.
- cols::Vector{Int64} Array containing the column indexes.
- Wdata::Vector{ComplexF64} Array containing the data for W, corresponding to the rows and cols.
- Idata::Vector{ComplexF64} Array containing the data for I, corresponding to the rows and cols.
- σ::Float64 Frequency for 'shift and invert', typically this is the TAE frequency.
- nev::Int64 Number of eigenvalues to solve for.
- n::Int64 Size of the matrix.
- dir::String Directory the results are written to.
"""
function old_par_solve(rows::Vector{Int64}, cols::Vector{Int64}, Wdata::Vector{ComplexF64}, Idata::Vector{ComplexF64}; σ::Float64, nev=5::Int64, n::Int64, dir::String)


    #eigenvalues are written in matlab format as this offers greater precision than default, 
    #and is easier to read from file
    evals_str = " -eps_view_values :" * dir * "vals.dat:ascii_matlab"
    #eigenfunctions are written as symmodu, which ignores which part of the eigenfunctions is stored by each proc.
    efuncs_str = " -eps_view_vectors :" * dir * "funcs.dat:ascii_symmodu"

    #these are combined with nev and σ, specifiying the setup for slepc.
    #eps_harmonics means we are searching inside the spectrum.
    #-st_type sinvert
    #slepcargs = @sprintf("-eps_nev %d -eps_target %s -eps_harmonic", nev, σ) * evals_str * efuncs_str
    slepcargs = @sprintf("-eps_nev %d -eps_target %s -st_type sinvert", nev, σ) * evals_str * efuncs_str
    

    SlepcInitialize(slepcargs)

    #The two empty matrices are created.
    W = create_matrix(n, n, auto_setup=true)
    I = create_matrix(n, n, auto_setup=true)

    #Values are added based on Coo format.
    set_values!(W, rows, cols, Wdata)
    set_values!(I, rows, cols, Idata)
    
    #matrices are then assembled.
    assemble!(W)
    assemble!(I)

    #we create the eps object, auto setup uses the slepcargs above.
    eps = create_eps(W, I; auto_setup=true)

    
    #solve the problem.
    #slepc args above will automatically write the solution to file.
    #note sometimes nev is ignored and solution just gives us the number of converged eigenvalues.
    solve!(eps)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        display(EPSGetConverged(eps))
    end
    
    #free the memory used.
    destroy!(W)
    destroy!(I)
    destroy!(eps)

    SlepcFinalize()
    
end
