"""
````
Grid=ExtendableGrids.simplexgrid
````
Re-Export of ExtendableGrids.simplexgrid
"""
const Grid = ExtendableGrids.simplexgrid


num_cellregions(grid::ExtendableGrid)=grid[NumCellRegions]
num_bfaceregions(grid::ExtendableGrid)=grid[NumBFaceRegions]

cellregions(grid::ExtendableGrid)= grid[CellRegions]
bfaceregions(grid::ExtendableGrid)= grid[BFaceRegions]
cellnodes(grid::ExtendableGrid)= grid[CellNodes]
bfacenodes(grid::ExtendableGrid)= grid[BFaceNodes]
coordinates(grid::ExtendableGrid)= grid[Coordinates]



function cartesian!(grid::ExtendableGrid)
    if dim_space(grid)==1
        grid[CoordinateSystem]=Cartesian1D
    elseif dim_space(grid)==2
        grid[CoordinateSystem]=Cartesian2D
    else dim_space(grid)==3
        grid[CoordinateSystem]=Cartesian3D
    end
    return grid
end


function circular_symmetric!(grid::ExtendableGrid)
    if dim_space(grid)==1
        grid[CoordinateSystem]=Polar1D
    elseif dim_space(grid)==2
        grid[CoordinateSystem]=Cylindrical2D
    else
        throw(DomainError(3,"Unable to handle circular symmetry for 3D grid"))
    end
    return grid
end

function spherical_symmetric!(grid::ExtendableGrid)
    d=dim_space(grid)
    if d==1
        grid[CoordinateSystem]=Spherical1D
    else
        throw(DomainError(d,"Unable to handle spherical symmetry for $(d)D grid"))
    end
    return grid
end


