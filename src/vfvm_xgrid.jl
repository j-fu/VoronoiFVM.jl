Grid(X)=simplexgrid(X)
Grid(X,Y)=simplexgrid(X,Y)
function Grid(;flags::String="pAaqDQ",
              points=Array{Cdouble,2}(undef,0,0),
              bfaces=Array{Cint,2}(undef,0,0),
              bfaceregions=Array{Cint,1}(undef,0),
              regionpoints=Array{Cdouble,2}(undef,0,0),
              regionnumbers=Array{Cint,1}(undef,0),
              regionvolumes=Array{Cdouble,1}(undef,0)
              )
    simplexgrid(flags=flags,
                points=points,
                bfaces=bfaces,
                bfaceregions=bfaceregions,
                regionpoints=regionpoints,
                regionnumbers=regionnumbers,
                regionvolumes=regionvolumes)
end



cellfactors!(grid::ExtendableGrid,icell,nodefac,edgefac)=cellfactors!(grid[CellGeometries][1],grid[CoordinateSystem],grid[Coordinates],grid[CellNodes], icell, nodefac,edgefac)
bfacefactors!(grid::ExtendableGrid,icell,nodefac)=bfacefactors!(grid[BFaceGeometries][1],grid[CoordinateSystem],grid[Coordinates],grid[CellNodes],icell,nodefac)



num_cellregions(grid::ExtendableGrid)=grid[NumCellRegions]
num_bfaceregions(grid::ExtendableGrid)=grid[NumBFaceRegions]


cellregions(grid::ExtendableGrid)= grid[CellRegions]
bfaceregions(grid::ExtendableGrid)= grid[BFaceRegions]
cellnodes(grid::ExtendableGrid)= grid[CellNodes]
bfacenodes(grid::ExtendableGrid)= grid[BFaceNodes]
coordinates(grid::ExtendableGrid)= grid[Coordinates]

num_nodes_per_cell(grid::ExtendableGrid)=num_targets(grid[CellNodes],1)
num_edges(grid::ExtendableGrid)=0 # !!!

num_nodes_per_bface(grid::ExtendableGrid)=num_targets(grid[CellNodes],1)-1
function num_edges_per_cell(grid::ExtendableGrid)
    nedges=[1,3,6]
    nedges[dim_space(grid)]
end

                   
reg_cell(grid::ExtendableGrid,i)=grid[CellRegions][i]
reg_bface(grid::ExtendableGrid,i)=grid[BFaceRegions][i]

cellnode(grid::ExtendableGrid,i,j)=grid[CellNodes][i,j]
bfacenode(grid::ExtendableGrid,i,j)=grid[BFaceNodes][i,j]

function celledgenode(grid::ExtendableGrid,inode,iedge,icell)
    cn=grid[CellNodes]

    if size(cn,1)==2
        local_celledgenodes=reshape([1 2],:,1)
    else
        # see grid/simplex.h in pdelib
        local_celledgenodes=zeros(Int64,2,3)
        local_celledgenodes[1,1]=2
        local_celledgenodes[2,1]=3

        local_celledgenodes[1,2]=3
        local_celledgenodes[2,2]=1
        
        local_celledgenodes[1,3]=1
        local_celledgenodes[2,3]=2
    end
    cn[local_celledgenodes[inode,iedge],icell]
end




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
