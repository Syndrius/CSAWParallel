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
    rind, θind, ζind, _ = index_to_grid(indstart+1, grids)

    grid_points = Array{Tuple{Int64, Int64, Int64}}(undef, npoints)
    n = 1
    #iterate through the grid points until all are acounted for.
    #when ζ reaches its max, θ is incremented
    #when θ reaches its max, r is incremented
    while n <= npoints
        grid_points[n] = (rind, θind, ζind)
        n += 1
        ζind += 1
        if ζind > grids.ζ.N
            ζind = 1
            θind += 1
            if θind > grids.θ.N
                θind = 1
                rind += 1
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

    #TODO 
    #this has been changed to match fff case, has not been verified yet!

    #note int32 is due to petsc.

    #total number of gridpoints controlled by this proc
    npoints = (indend - indstart) ÷ 4 #4 for ffs
    
    #gets the first grid point controlled by this proc
    #+1 to convert to julia 
    rind, θind, _, _ = index_to_grid(indstart+1, grids)

    grid_points = Array{Tuple{Int64, Int64}}(undef, npoints)
    n = 1
    #for this case, the grid is only divided by r, θ, all ζ values will always be on the same proc for a given r, θ
    #iterate through the grid points until all are acounted for.
    #when θ reaches its max, r is incremented
    while n <= npoints
        grid_points[n] = (rind, θind)
        n += 1
        θind += 1
        if θind > grids.θ.N
            θind = 1
            rind += 1
        end
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

