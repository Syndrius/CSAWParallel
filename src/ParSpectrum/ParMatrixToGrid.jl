
#converting the matrix ownership range into the actual grids we want to loop over
#need to do this for the other cases as well.
function matrix_to_grid(indstart::Int32, indend::Int32, grids::MID.FFFGridsT)

    
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
    
    return grid_points


end