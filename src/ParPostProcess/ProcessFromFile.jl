
"""
    par_post_process(dir::String)

Post processes the parallel solutions that have been written to file.
This is done outisde the main functions as there is a bug when converting petsc vec's to julia arrays.
This function is run with a single proc.
"""
function par_post_process(dir::String, deriv=false)
    mkpath(joinpath(dir, "efuncs"))
    mkpath(joinpath(dir, "efuncs_ft"))

    prob, grids, _ = inputs_from_file(dir=dir)

    vals = load_object(joinpath(dir, "vals_raw.jld2"))
    nevals = length(vals)

    errs = load_object(joinpath(dir, "errs.jld2"))

    #tolerance used by slepc.
    tol = errs[end]
    errs = errs[1:end-1]

    #creates arrays to store each solution while being processed
    ϕp, ϕpft = PostProcessing.allocate_phi_arrays(grids, deriv=deriv)

    #creates the appropriate fourier transform plan.
    plan = PostProcessing.create_ft_plan(ϕpft, grids)

    x1grid = inst_grids(grids)[1]

    #arrays to store the maximum value of ϕft and the correspondning r value.
    rmarray = Array{Int64}(undef, grids.x2.N, grids.x3.N)
    ϕmarray = Array{Float64}(undef, grids.x2.N, grids.x3.N)

    
    #stores the unique inds
    un_inds = Int64[]

    #first we determine the unique evals.
    #most relevant for slice solving where we might overlap the spectrum.
    for i in 1:nevals

        unique_found = true #assume this one is unique.

        for j in un_inds
            #eigenvalue matches previosuly found eigenvalue.
            if abs(vals[j] - vals[i]) < 10*tol * abs(vals[i])

                #check to see if the eigenfunctions are approximatly the same
                if equal_efuncs(dir, i, j)
                    if errs[j] > errs[i]
                        #the jth index is replaced by i
                        #as this eigenfunctions is the same but with lower error.
                        un_inds[un_inds .== j] .= i
                        unique_found = false
                        #otherwise the jth, i.e. already stored ind is better
                    else
                        unique_found = false
                    end
                end
            end
        end
        if unique_found
            push!(un_inds, i)
        end
    end

    un_nevals = length(un_inds)
    x1ms = Array{Float64}(undef, un_nevals)
    ω = zeros(ComplexF64, un_nevals)
    mode_labs = Tuple{Int, Int}[] 

    #unique eigenvalues are processed and written to file.
    for i in 1:un_nevals

        efunc_read = @sprintf("efunc%05d.hdf5", un_inds[i]) 
        efunc_write = @sprintf("efunc%05d.jld2", i)

        #unfort doesn't handle complex numbers v well
        efunc_split = load_object(joinpath(dir, "efuncs_raw", efunc_read))

        #so we manually create the solution properly
        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im

        #solution is reconstructed into the 3d grid.
        PostProcessing.reconstruct_phi!(efunc, grids, ϕp, ϕpft, plan)
        
        #solution is labeled by the largest harmonic
        x1ind, mode_lab = PostProcessing.label_mode(ϕpft, grids, rmarray, ϕmarray)

        push!(mode_labs, mode_lab)

        #storing the radial location of the peak, for the continuum
        x1ms[i] = x1grid[x1ind]

        #normalise the eigenvalues
        ω[i] = prob.geo.R0 * sqrt(vals[un_inds[i]])

        save_object(joinpath(dir, "efuncs", efunc_write), ϕp)
        save_object(joinpath(dir, "efuncs_ft", efunc_write), ϕpft)
        
    end
    evals = EvalsT(ω, x1ms, mode_labs)

    save_object(joinpath(dir, "evals.jld2"), evals)
    save_object(joinpath(dir, "unique_inds.jld2"), un_inds)
end


"""
    equal_efuncs(dir::String, ind1::Int64, ind2::Int64)

Checks if two eigenfunctions are equivalent.
"""
function equal_efuncs(dir::String, ind1::Int64, ind2::Int64)
    efunc_read1 = @sprintf("efunc%05d.hdf5", ind1)
    
    efunc_split = load_object(joinpath(dir, "efuncs_raw", efunc_read1))
    efunc1 = efunc_split[1, :] .+ efunc_split[2, :] * 1im

    efunc_read2 = @sprintf("efunc%05d.hdf5", ind2)
    efunc_split = load_object(joinpath(dir, "efuncs_raw", efunc_read2))
    efunc2 = efunc_split[1, :] .+ efunc_split[2, :] * 1im
    #computes the inner product of the two and each individually to see if they are the same.
    if 0.99 < abs(dot(efunc1, efunc2))^2 / abs(dot(efunc1, efunc1)) / abs(dot(efunc2, efunc2))
        return true
    end
    return false
end

