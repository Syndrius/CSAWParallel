
#converting the matrix ownership range into the actual grids we want to loop over
#need to do this for the other cases as well.
function matrix_to_grid(indstart::Int64, indend::Int32, grids::MID.FFFGridsT)

    
    #this is going to be v tru=icky, as θstart will be different for rstart
    #otherwise θstart will just be 1. Same for ζ, although ζ split is unlikely I think.
    #May need to make main loop into a function so that it can be called multiple times
    #will be a bit annoying though. The function will need like a bajillion args...
    #other options is we consider all tuples of r, θ, ζ
    #think that is a better option.
    rstart = div(indstart-1, 8*grids.θ.N*grids.ζ.N) + 1
    θstart = mod(div(indstart-1, 8*grids.ζ.N), grids.θ.N) + 1
    ζstart = mod(div(indstart-1, 8), grids.ζ.N) + 1

    rend = div(indend-1, 8*grids.θ.N*grids.ζ.N) + 1
    θend = mod(div(indend-1, 8*grids.ζ.N), grids.θ.N) + 1
    ζend = mod(div(indend-1, 8), grids.ζ.N) + 1


    @printf("Indstart of %d, gives (%d, %d, %d)\n", indstart, rstart, θstart, ζstart)
    @printf("Indend of %d, gives (%d, %d, %d)\n", indend, rend, θend, ζend)

    grid_points = Tuple{Int, Int, Int}[]

    for i in ζstart:grids.ζ.N
        push!(grid_points, (rstart, θstart, i))
    end

    for i in θstart+1:grids.θ.N

        for j in 1:grids.ζ.N
            push!(grid_points, (rstart, i, j))
        end
    end

    for i in rstart+1:rend-1
        for j in 1:grids.θ.N
            for k in 1:grids.ζ.N
                push!(grid_points, (i, j, k))
            end
        end
    end

    for i in 1:θend-1
        for j in 1:grids.ζ.N
            push!(grid_points, (rend, i, j))
        end
    end

    for i in 1:ζend
        push!(grid_points, (rend, θend, i))
    end

    @printf("For indstart, indend of (%d, %d), we have (%d)\n", indstart, indend, size(grid_points)[1])
    
    #display(rstart)
    #display(θstart)
    #display(ζstart)
    #display(rend)
    #display(θend)
    #display(ζend)
    #display(grid_points)
    return grid_points


end

#TODO
#this,
#the first int is 64 because it is modified, perhaps we should modify it in here instead???
function matrix_to_grid(indstart::Int64, indend::Int32, grids::MID.FFSGridsT)

    #for this case, the grid is only divided by r, θ, all ζ values will always be on the same proc for a given r, θ
    rstart, θstart, _, _ = index_to_grid(indstart, grids)
    rend, θend, _, _ = index_to_grid(Int64(indend), grids)

    grid_points = Tuple{Int, Int}[]

    #handles cases where where proc only has a portion of the θgrid for rstart
    for i in θstart:grids.θ.N
        push!(grid_points, (rstart, i))
    end

    #normal casees, which have all of r and all or θ
    for i in rstart+1:rend-1
        for j in 1:grids.θ.N
            push!(grid_points, (i, j))
        end
    end

    #handles cases where proc only has portion of θgrid for rend
    for i in 1:θend
        push!(grid_points, (rend, i))
    end

    return grid_points

end


function matrix_to_grid(indstart::Int32, indend::Int32, grids::MID.FSSGridsT)

    #this doesn't really make sense for this case!

    #this should be v different.
    #ideally we will still have the same structure, ie loop
    #through grid_points.

    rstart, θstart, ζstart, _ = MID.index_to_grid(indstart, grids)
    rend, θend, ζend, _ = MID.index_to_grid(indend, grids)


    grid_points = Tuple{Int, Int, Int}[]

    for i in ζstart:grids.ζ.count
        push!(grid_points, (rstart, θstart, i))
    end

    for i in θstart+1:grids.θ.count

        for j in 1:grids.ζ.N
            push!(grid_points, (rstart, i, j))
        end
    end

    for i in rstart+1:rend-1
        for j in 1:grids.θ.count
            for k in 1:grids.ζ.count
                push!(grid_points, (i, j, k))
            end
        end
    end

    for i in 1:θend-1
        for j in 1:grids.ζ.count
            push!(grid_points, (rend, i, j))
        end
    end

    for i in 1:ζend
        push!(grid_points, (rend, θend, i))
    end
    
    return grid_points


end
