"""
    matrix_to_grid(indstart::Int32, indend::Int32, grids::FFFGridsT)

Converts the ownership range of the matrix into grid points for the core to iterate over.
"""
function matrix_to_grid(indstart::Int32, indend::Int32, grids::FFFGridsT)
    
    #note int32 is due to petsc.

    #total number of gridpoints controlled by this proc
    npoints = (indend - indstart) ÷ 8 #8 for fff
    
    #gets the first grid point controlled by this proc
    #+1 to convert to julia
    x1ind, x2ind, x3ind, _ = index_to_grid(indstart+1, grids)

    grid_points = Array{Tuple{Int64, Int64, Int64}}(undef, npoints)
    n = 1
    #iterate through the grid points until all are acounted for.
    #when x3 reaches its max, x2 is incremented
    #when x2 reaches its max, r is incremented
    while n <= npoints
        grid_points[n] = (x1ind, x2ind, x3ind)
        n += 1
        x3ind += 1
        if x3ind > grids.x3.N
            x3ind = 1
            x2ind += 1
            if x2ind > grids.x2.N
                x2ind = 1
                x1ind += 1
            end
        end
    end
    return grid_points

end


"""
    matrix_to_grid(indstart::Int32, indend::Int32, grids::FFFGridsT)

Converts the ownership range of the matrix into grid points for the core to iterate over.
"""
function matrix_to_grid(indstart::Int32, indend::Int32, grids::FFSGridsT)

    #note int32 is due to petsc.

    #total number of gridpoints controlled by this proc
    npoints = (indend - indstart) ÷ 4 #4 for ffs
    
    #gets the first grid point controlled by this proc
    #+1 to convert to julia 
    x1ind, x2ind, _, _ = index_to_grid(indstart+1, grids)

    grid_points = Array{Tuple{Int64, Int64}}(undef, npoints)
    n = 1
    #for this case, the grid is only divided by r, x2, all x3 values will always be on the same proc for a given r, x2
    #iterate through the grid points until all are acounted for.
    #when x2 reaches its max, r is incremented
    while n <= npoints
        grid_points[n] = (x1ind, x2ind)
        n += 1
        x2ind += 1
        if x2ind > grids.x2.N
            x2ind = 1
            x1ind += 1
        end
    end
    return grid_points
end


"""
    matrix_to_grid(indstart::Int32, indend::Int32, grids::FFFGridsT)

Converts the ownership range of the matrix into grid points for the core to iterate over.
"""
function matrix_to_grid(indstart::Int32, indend::Int32, grids::FSSGridsT)

    #for this case, the grid is only divided by r, all x2, x3 values will always be on the same proc for a given r
    #+1 to convert to julia 
    x1start, _, _, _ = index_to_grid(indstart+1, grids)
    #+1 to change to julia, but then -1 to use julia's inclusive indexing
    x1end, _, _, _ = index_to_grid(Int64(indend), grids)

    return collect(x1start:x1end)
end 

