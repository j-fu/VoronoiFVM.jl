#=
Definition of subgrid interface methods
=#

##########################################################
"""
$(TYPEDSIGNATURES)

Space dimension of grid
"""
dim_space(grid::SubGrid)= size(grid.coord,1)


##########################################################
"""
$(TYPEDSIGNATURES)


Number of nodes in grid
"""
num_nodes(grid::SubGrid)= size(grid.coord,2)


##########################################################
"""
$(TYPEDSIGNATURES)

Number of cells in grid
"""
num_cells(grid::SubGrid)= size(grid.cellnodes,2)

##########################################################
"""
$(TYPEDSIGNATURES)

Number of edges in grid
"""
num_edges(grid::SubGrid)= size(grid.edgenodes,2)

##########################################################
"""
$(TYPEDSIGNATURES)

Return index of i-th local node in cell icell
"""
cellnode(grid::SubGrid,inode,icell)=grid.cellnodes[inode,icell]

##########################################################
"""
$(TYPEDSIGNATURES)

Return view of coordinates of node `inode`.
"""
nodecoord(grid::SubGrid,inode)=view(grid.coord,:,inode)

##########################################################
"""
$(TYPEDSIGNATURES)

Return number of nodes per cell in grid.
"""
num_nodes_per_cell(grid::SubGrid)= size(grid.cellnodes,1)

##########################################################
"""
$(TYPEDSIGNATURES)

Return element type of grid coordinates.
"""
Base.eltype(grid::SubGrid)=Base.eltype(grid.coord)

