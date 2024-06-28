""" 

Handles Input and Output for parallel cases. Most notable is reading the solutions from file, which have weird format.


"""
module ParIo

using MID
using MIDParallel.ParSpectrum

using DelimitedFiles #not certain this will be used!
using Printf

export par_construct_to_file
export par_solve_from_file
export par_spectrum_from_file

export par_vals_from_file
export par_funcs_from_file
export par_func_from_file


"""
Computes the spectrum in parallel based on inputs read from file.

# Args
- dir::String Directory the inputs are stored and results are written.
- freq::Float64 Frequency for 'shift and invert', typically this is the TAE frequency.
- nev::Int64 Number of eigenvalues to solve for.
"""
function par_spectrum_from_file(; dir::String, freq::Float64, nev=5::Int64)

    prob, grids = inputs_from_file(dir=dir)

    par_construct_and_solve(prob=prob, grids=grids, σ=freq, nev=nev, dir=dir)

end


#this is probably a useless function now!
#could be modified in future if we ever want to solve different parts of the spectrum?
#=
function par_solve_from_file(; dir::String)


    file_inds = dir * @sprintf("inds.dat")
    file_data = dir * @sprintf("data.dat")

    rows, cols = eachrow(readdlm(file_inds, ',', Int64))
    #diffucult to have access to n when solving from file, test just doing this for now
    #this will hopefully change when the problem is written to file etc.
    n = maximum(cols) #this is giving something cooked???? Could be a srs worry.
    n=400 #hard code for now, needs to change
    R0 = 10.0
    tae_freq = (0.381 / R0)^2
    Wdata, Idata = eachrow(readdlm(file_data, ',', ComplexF64))
    MPI.Init()
    #this is almost certainly solving the same matrix on each proc.
    par_solve(rows, cols, Wdata, Idata, R0=R0, n=n, dir=dir, σ=tae_freq)

    MPI.Finalize() 

end
=#

"""
Constructs the matrix in parallel then writes the matrix to file.
Most likley to be used when parallel solver is not functional.

# Args
- prob::ProblemT Struct containing the functions and parameters that define the problem we are solving
- grids::GridT Grids to solve over.
dir::String Directory for matrices to be written to.
"""
function par_construct_to_file(; prob::MID.ProblemT, grids::MID.GridsT, dir::String)

    MPI.Init()

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    root = 0
    rows, cols, Wdata, Idata = par_construct(prob=prob, grids=grids)
    MPI.Barrier(comm) #not sure if this is needed!
    #for some reason MPI used Int32's
    global_counts = Int32.(MPI.Allgather([length(rows)], comm))

    global_rows = MPI.Gatherv(rows, global_counts, root, comm)
    global_cols = MPI.Gatherv(cols, global_counts, root, comm)
    global_Wdata = MPI.Gatherv(Wdata, global_counts, root, comm)
    global_Idata = MPI.Gatherv(Idata, global_counts, root, comm)
    
    #I think this would be better if it was written as a single file with 4 columns
    if rank == root
        #may want to make the file more sophisticated so that we separate the file from the solutions.
        file_name = dir*@sprintf("inds.dat")
        open(file_name, "w") do file
            writedlm(file, [global_rows, global_cols], ",")
        end 

        file_name = dir*@sprintf("data.dat")
        open(file_name, "w") do file
            writedlm(file, [global_Wdata, global_Idata], ",")
        end 
    end
    MPI.Finalize()

end

"""
Reads the eigenvalues from file after they have been output by the parallel solver.

# Args
- filname::String File where eigenvalues are stored.
- R0::Float64 Major Radius, for nomalising.
"""
function par_vals_from_file(filename::String, R0::Float64)
    #efuncs are written to 'matlab' formatted files by slepc. 
    #This includes some header information that needs ot be removed.
    evals = open(filename, "r") do file
        s = readlines(file)[2:end-1] 
        parse.(ComplexF64, s)
    end
    #returns the normalised form.
    return R0 .* sqrt.(evals)
end


"""
Reads the eigenfunctions from file after they have been output by the parallel solver. Returns the reconstructed form of the eigenfunctions.

# Args
- filname::String File where eigenfunctions are stored.
- nevals::Int64 Number of eigenvalues and eigenvectors.
- nr::Int64 Number of radial grid points, used for reconstructing.
- nm::Int64 Number of poloidal modes, used for reconstructing.
- nn::Int64 Number of toroidal modes, used for reconstructing.
"""
function par_funcs_from_file(filename::String, nevals::Int64, grids::GridsT)

    n = matrix_dim(grids) #size of matrix, note reconstruct phi removes the derivatives.

    efuncs = open(filename, "r") do file
        efuncs = Array{ComplexF64}(undef, n, nevals) 
        #temp array for combining two floats into a single complex number.
        tmp = Array{ComplexF64}(undef, n) 
        s = readlines(file)
        for i in 1:nevals
            
            #awkward pattern avoids Slepc header before each vector.
            for (j, str) in enumerate(split.(s[3+(i-1)*(n+2):2*i+i*n]))
                real = parse(Float64, str[1])
                imag = parse(Float64, str[2]) #* 1im
                tmp[j] = real + imag*1im
            end
            efuncs[:, i] = tmp
        end
    
        return efuncs
    end

    return reconstruct_phi(efuncs, nevals, grids)

end



"""
Reads a single eigenfunction from file after they have been output by the parallel solver. Returns the reconstructed form of the eigenfunction. Much faster and more memory efficient than reading all eigenfucntions.

# Args
- filname::String File where eigenfunctions are stored.
- ind::Int64 Index of desired eigenfunction.
- nr::Int64 Number of radial grid points, used for reconstructing.
- nm::Int64 Number of poloidal modes, used for reconstructing.
- nn::Int64 Number of toroidal modes, used for reconstructing.
"""
function par_func_from_file(filename::String, ind::Int64, grids::GridsT)

    n = matrix_dim(grids) #size of matrix, note reconstruct phi removes the derivatives.

    efunc = open(filename, "r") do file
        efunc = Array{ComplexF64}(undef, n, 1) 

        s = readlines(file)
        for (i, str) in enumerate(split.(s[3+(ind-1)*(n+2):2*ind+ind*n]))
            real = parse(Float64, str[1])
            imag = parse(Float64, str[2]) #* 1im
            efunc[i, 1] = real + imag*1im
        end
        
        return efunc
    end

    return reconstruct_phi(efunc, 1, grids)

end


end