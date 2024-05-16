
module Parallel

#so this works, obvs without registering MID, we have to add this locally, ie Pkg.add ~/phd/MID or whatever
using MID
using MPI
using FFTW
using FastGaussQuadrature
using SparseArrays
#using LinearAlgebra
#using Arpack
using DelimitedFiles #not certain this will be used!
using Printf


#note to get these files to work on mac, had to modify the load.jl files in both cases
#mac uses dylib files, while linux uses .so files so PetscWrap and SlepcWrap were unable to find the files.
using PetscWrap
using SlepcWrap



include("ParConstruct.jl")

export par_construct


include("ParSolve.jl")

export par_solve

export par_construct_and_solve
export par_construct_to_file
export par_solve_from_file


#cannot get this to work, Julia MPI defaults to MPICH, and cannot easily be changed, Petsc is unable to install with MPICH for some reason.
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

#for writing shit straight to file so we can avoid many of the issues.
function par_construct_to_file(; prob::MID.ProblemT, grids::MID.GridsT, dir::String)

    MPI.Init()

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    root = 0
    #something may be cooked here!
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


#doesn't allow reconstruction, as that will be v difficult!
#may want option to write the matrix to file as well? Not certain it will be that useful, but maybe being able to solve in different cases will be handy?
function par_construct_and_solve(; prob::MID.ProblemT, grids::MID.GridsT, efuncs=true::Bool, σ=0.0::Float64, nev=20::Int64, dir::String)

    MPI.Init()

    root = 0

    comm = MPI.COMM_WORLD

    #this does not work atm as, with n=2 for example, rather than consructing the global matrix,
        #it constructs some stupid af combo of 2 smaller matrices.
        # so it seems that the construct thingo needs to know what all the data is. 
        
    rows, cols, Wdata, Idata = par_construct(prob=prob, grids=grids)

    global_counts = Int32.(MPI.Allgather([length(rows)], comm))

    #this seems like a sub optimal solution, given this requires all procs have the same memory???
    #But this does work, and is essentially the same as reading from file, but we skip the read-write step
    global_rows = MPI.Allgatherv(rows, global_counts, comm)
    global_cols = MPI.Allgatherv(cols, global_counts, comm)
    global_Wdata = MPI.Allgatherv(Wdata, global_counts, comm)
    global_Idata = MPI.Allgatherv(Idata, global_counts, comm)

    n = 2 * grids.rd.N * grids.pmd.count * grids.tmd.count
    par_solve(global_rows, global_cols, global_Wdata, global_Idata, σ=σ, nev=nev, R0=prob.geo.R0, n=n, dir=dir)

    MPI.Finalize()

end

#bad name!
function par_spectrum_from_file(; dir::String, freq::Float64, nev=20::Int64)

    prob, grids = inputs_from_file(dir=dir)

    par_construct_and_solve(prob=prob, grids=grids, σ=freq, nev=nev, dir=dir)

end

end