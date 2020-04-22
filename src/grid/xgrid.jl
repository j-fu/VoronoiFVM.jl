coord_type(grid::ExtendableGrid)=Base.eltype(grid[Coordinates])
index_type(grid::ExtendableGrid)=Base.eltype(grid[CellNodes])

num_cellregions(grid::ExtendableGrid)=grid[NumCellRegions]
num_bfaceregions(grid::ExtendableGrid)=grid[NumBFaceRegions]
num_nodes(grid::ExtendableGrid)= size(grid[Coordinates],2)
num_cells(grid::ExtendableGrid)= num_sources(grid[CellNodes])
num_bfaces(grid::ExtendableGrid)= num_sources(grid[BFaceNodes])
cellregions(grid::ExtendableGrid)= grid[CellRegions]
bfaceregions(grid::ExtendableGrid)= grid[BFaceRegions]
cellnodes(grid::ExtendableGrid)= grid[CellNodes]
bfacenodes(grid::ExtendableGrid)= grid[BFaceNodes]
dim_grid(grid::ExtendableGrid)= dim_space(grid)-1
dim_space(grid::ExtendableGrid)= size(grid[Coordinates],1)
num_nodes_per_cell(grid::ExtendableGrid)=num_targets(grid[CellNodes],1)
num_edges(grid::ExtendableGrid)=0 # !!!

num_nodes_per_bface(grid::ExtendableGrid)=num_targets(grid[CellNodes],1)-1
function num_edges_per_cell(grid::ExtendableGrid)
    nedges=[1,3,6]
    nedges[dim_space(grid)]
end

cellfactors!(grid::ExtendableGrid,icell,nodefac,edgefac)=cellfactors!(grid[CellTypes][1],grid[CoordinateSystem],grid[Coordinates],grid[CellNodes], icell, nodefac,edgefac)
bfacefactors!(grid::ExtendableGrid,icell,nodefac)=bfacefactors!(grid[BFaceTypes][1],grid[CoordinateSystem],grid[Coordinates],grid[CellNodes],icell,nodefac)
                   
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
        local_celledgenodes=zeros(Ti,2,3)
        local_celledgenodes[1,1]=2
        local_celledgenodes[2,1]=3

        local_celledgenodes[1,2]=3
        local_celledgenodes[2,2]=1
        
        local_celledgenodes[1,3]=1
        local_celledgenodes[2,3]=2
    end
    cn[local_celledgenodes[inode,iedge],icell]
end
