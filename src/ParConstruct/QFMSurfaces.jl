
#this file computes the qfm surfaces in parallel, simply call the QFM Surface function, but the list of rational surfaces is split across processors

#perhaps the name of this should be different!

#TODO, hasn't been updated to new surface constructtion

"""

Computes the qfm surfaces in parallel by divided them up between cores.
Currently there is an issue broadcasting the surface array, so they are written to file individually.
"""
function par_construct_surfaces(plist::Array{Int64}, qlist::Array{Int64}, sguesslist::Array{Float64}, prob::ProblemT, dir::String, MM=4::Int64, M=24::Int64, N=8::Int64)


    #think this function is designed to be called in its own session
    MPI.Init()

    root = 0 

    met = MetT()
    B = BFieldT()

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    nprocs = MPI.Comm_size(comm) #total number of workers including root.

    local_plist = mpi_split_array(plist)
    local_qlist = mpi_split_array(qlist)
    local_slist = mpi_split_array(sguesslist)


    local_surfs = construct_surfaces(local_plist, local_qlist, local_slist, prob)

    #so the isbitstype is causing lots of issues, just going to write each proc individually, then they can be combined later
    save_object(dir * "local_surfs" * string(rank) * ".jld2", local_surfs)

    MPI.Finalize()
    
end



"""
    mpi_split_array(array::AbstractArray)

Function to split an array evenly between cores.
"""
function mpi_split_array(array::AbstractArray)
    #assumes MPI.init() has already been called

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm) #total number of workers including root.
    
    counts = zeros(Int64, nprocs)
    splits = zeros(Int64, nprocs+1)
    ar_size = length(array)

    counts_guess = Int64(div(ar_size, nprocs, RoundDown))
    Remainder = Int64(ar_size - counts_guess*nprocs)
    counts .= counts_guess
    for i in 1:Remainder
        counts[i] += 1
    end

    splits[1:end-1] .= cumsum(append!([0], counts))[1:nprocs] .+ 1
    splits[end] = ar_size+1


    return array[splits[rank+1]:(splits[rank+2] -1)]
end


"""
    gather_surfs(dir::String, procs::Int64)

Temporary function until we can share the surf struct with MPI.
This function just gathers the surfaces written to file by each core into a single struture.
"""
function gather_surfs(dir::String, procs::Int64)

    surfs = QFMSurfaceT[]

    for i in 1:procs
        local_surfs = load_object(dir * "local_surfs"*string(i-1) * ".jld2")

        for j in local_surfs
            push!(surfs, j)
        end
    end
    save_object(dir * "qfm_benchmark_surfaces.jld2", surfs)
end
