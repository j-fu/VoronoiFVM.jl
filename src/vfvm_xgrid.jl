"""
````
Grid=ExtendableGrids.simplexgrid
````
Re-Export of ExtendableGrids.simplexgrid
!!! compat 
    Deprecated as of v1.20
"""
const Grid = ExtendableGrids.simplexgrid

@deprecate Grid ExtendableGrids.simplexgrid

cellregions(grid::ExtendableGrid) = grid[CellRegions]
bfaceregions(grid::ExtendableGrid) = grid[BFaceRegions]
cellnodes(grid::ExtendableGrid) = grid[CellNodes]
bfacenodes(grid::ExtendableGrid) = grid[BFaceNodes]

"""
    coordinates(grid::ExtendableGrid)

Return coordinate array of grid. 

!!! compat 
    Deprecated as of v1.20
"""
coordinates(grid::ExtendableGrid) = grid[Coordinates]

@deprecate coordinates(grid) grid[Coordinates]

"""
   $(SIGNATURES)

Set coordinate system in grid to cartesian.
"""
function cartesian!(grid::ExtendableGrid)
    if dim_space(grid) == 1
        grid[CoordinateSystem] = Cartesian1D
    elseif dim_space(grid) == 2
        grid[CoordinateSystem] = Cartesian2D
    else
        dim_space(grid) == 3
        grid[CoordinateSystem] = Cartesian3D
    end
    return grid
end

"""
   $(SIGNATURES)

Set coordinate system in grid to circular/cylindrical symmetry.
"""
function circular_symmetric!(grid::ExtendableGrid)
    if dim_space(grid) == 1
        grid[CoordinateSystem] = Polar1D
    elseif dim_space(grid) == 2
        grid[CoordinateSystem] = Cylindrical2D
    else
        throw(DomainError(3, "Unable to handle circular symmetry for 3D grid"))
    end
    return grid
end

"""
   $(SIGNATURES)

Set coordinate system in grid to spherical symmetry.
"""
function spherical_symmetric!(grid::ExtendableGrid)
    d = dim_space(grid)
    if d == 1
        grid[CoordinateSystem] = Spherical1D
    else
        throw(DomainError(d, "Unable to handle spherical symmetry for $(d)D grid"))
    end
    return grid
end
