"""
    matrix_to_grid(indstart::Int32, indend::Int32, grids::FFFGridsT)

Converts the ownership range of the matrix into grid points for the core to iterate over.
"""
function matrix_to_grid(indstart::Int32, indend::Int32, grids::FFFGridsT)
    
    #note int32 is due to petsc.

    #gets the start and final points for each coordinate.

    #+1 to convert to julia
    rstart, θstart, ζstart, _ = index_to_grid(indstart+1, grids)
    #+1 to change to julia, but then -1 to use julia's inclusive indexing
    rend, θend, ζend, _ = index_to_grid(Int64(indend), grids)


    grid_points = Tuple{Int64, Int64, Int64}[]

    #we now have to consider the edge of each block.
    #i.e. the grid may be split so that for the first θ point only some of the ζrange 
    #is included for this core.

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


"""
    matrix_to_grid(indstart::Int32, indend::Int32, grids::FFFGridsT)

Converts the ownership range of the matrix into grid points for the core to iterate over.
"""
function matrix_to_grid(indstart::Int32, indend::Int32, grids::FFSGridsT)

    #for this case, the grid is only divided by r, θ, all ζ values will always be on the same proc for a given r, θ
    #+1 to convert to julia 
    rstart, θstart, _, _ = index_to_grid(indstart+1, grids)
    #+1 to change to julia, but then -1 to use julia's inclusive indexing
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


"""
    matrix_to_grid(indstart::Int32, indend::Int32, grids::FFFGridsT)

Converts the ownership range of the matrix into grid points for the core to iterate over.
"""
function matrix_to_grid(indstart::Int32, indend::Int32, grids::FSSGridsT)

    #for this case, the grid is only divided by r, all θ, ζ values will always be on the same proc for a given r
    #+1 to convert to julia 
    rstart, _, _, _ = index_to_grid(indstart+1, grids)
    #+1 to change to julia, but then -1 to use julia's inclusive indexing
    rend, _, _, _ = index_to_grid(Int64(indend), grids)

    return collect(rstart:rend)
end 

