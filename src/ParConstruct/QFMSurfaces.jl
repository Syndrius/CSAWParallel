
#this file computes the qfm surfaces in parallel, simply call the QFM Surface function, but the list of rational surfaces is split across processors

#perhaps the name of this should be different!


#ok so this now works, bit of a disaster, but does the job
function par_construct_surfaces(plist, qlist, sguesslist, prob, dir)

    #this is still very good
    MM = 4
    M = 24
    N = 8

    #think this function is designed to be called in its own session
    MPI.Init()

    root = 0 #or 1?

    met = MetT()
    B = BFieldT()

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)


    nprocs = MPI.Comm_size(comm) #total number of workers including root.
    #I cannt beleive this does this, fkn useless.
    #local_plist = collect(Iterators.partition(plist, nprocs))[rank+1]
    #local_qlist = collect(Iterators.partition(qlist, nprocs))[rank+1]
    #local_slist = collect(Iterators.partition(sguesslist, nprocs))[rank+1]

    #this is a suboptimal split, the last proc will be doing less work, but who cares tbh.
    #local_plist = collect(Iterators.partition(plist, length(plist) ÷ nprocs + 1))[rank+1]
    #local_qlist = collect(Iterators.partition(qlist, length(qlist) ÷ nprocs + 1))[rank+1]
    #local_slist = collect(Iterators.partition(sguesslist, length(sguesslist) ÷ nprocs + 1))[rank+1]
    local_plist = mpi_split_array(plist)
    local_qlist = mpi_split_array(qlist)
    local_slist = mpi_split_array(sguesslist)


    local_surfs = construct_surfaces(local_plist, local_qlist, local_slist, prob)

    #so we cant just send the surfaces because they are a weird type. Instead we will send all of the data individually!
    #surfaces = MPI.Gather(local_surfs, root, comm)
    #=
    nlocal_surfs = length(local_surfs)

    #you would think that we do actually know the size of these by now!
    #we should change the name of these now we now what they actually are!
    local_ρ = zeros(nlocal_surfs)
    local_scos = zeros(nlocal_surfs, M+1, 2 * N + 1)
    local_tsin = zeros(nlocal_surfs, M+1, 2 * N + 1)
    local_ssin = zeros(nlocal_surfs, M+1, 2 * N + 1)
    local_tcos = zeros(nlocal_surfs, M+1, 2 * N + 1)

    local_scos = Array{Array{Float64, 2}, 1}
    local_tsin = Array{Array{Float64, 2}, 1}
    local_ssin = Array{Array{Float64, 2}, 1}
    local_tcos = Array{Array{Float64, 2}, 1}

    local_scos = Array{Float64, 2}[]
    local_tsin = Array{Float64, 2}[]
    local_ssin = Array{Float64, 2}[]
    local_tcos = Array{Float64, 2}[]


    for i in 1:nlocal_surfs
        local_ρ[i] = local_surfs[i].ρ
        #display(local_surfs[i].scos)
        push!(local_scos, local_surfs[i].scos)
        push!(local_tsin, local_surfs[i].tsin)
        push!(local_ssin, local_surfs[i].ssin)
        push!(local_tcos, local_surfs[i].tcos)
    end


    global_ρ = MPI.Gather(local_ρ, root, comm)
    global_scos = MPI.Gather(local_scos, root, comm)
    global_tsin = MPI.Gather(local_tsin, root, comm)
    global_ssin = MPI.Gather(local_ssin, root, comm)
    global_tcos = MPI.Gather(local_tcos, root, comm)

    if rank == root
        global_surfs = MID.QFMSurfaceT[]
        for i in 1:length(global_ρ)
            push!(global_surfs, MID.QFMSurfaceT(global_ρ[i], global_scos[i], global_tsin[i], global_ssin[i], global_tcos[i]))
        end
        #perhaps this should add jld2 on the end?
        display(length(global_surfs))
        save_object(filename, global_surfs)
    end
    =#
    #so the isbitstype is causing lots of issues, just going to write each proc individually, then they can be combined later

    save_object(dir * "local_surfs" * string(rank) * ".jld2", local_surfs)

    MPI.Finalize()
    
end


function mpi_split_array(array)
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
    #display(counts)

    splits[1:end-1] .= cumsum(append!([0], counts))[1:nprocs] .+ 1
    splits[end] = ar_size+1

    #display(splits)

    return array[splits[rank+1]:(splits[rank+2] -1)]
end

#stupid function because we couldn't get MPI gather to work for the structs storing the qfm surfaces.
#this is a disaster, as we cannot call this easily, however, it does work!
function gather_surfs(dir, procs)

    surfs = QFMSurfaceT[]

    for i in 1:procs
        local_surfs = load_object(dir * "local_surfs"*string(i-1) * ".jld2")

        for j in local_surfs
            push!(surfs, j)
        end
    end
    save_object(dir * "qfm_benchmark_surfaces.jld2", surfs)
end
