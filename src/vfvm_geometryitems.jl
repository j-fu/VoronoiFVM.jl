abstract type AbstractGeometryItem end

##################################################################
"""
$(TYPEDEF)

Structure holding local boundary  node information.

$(TYPEDFIELDS)
"""
mutable struct BNode{Tv} <: AbstractGeometryItem 

    """
    Index in grid
    """
    index::Int32

    """
    Boundary region number
    """
    region::Int32

    """
    Number of species defined in node
    """
    nspec::Int64

    """
    Grid
    """
    grid
    
    """
    Physics data
    """
    physics_data

    """
    (deprecated) Coordinates
    """
    coord::Array{Tv,1}

    
    BNode{Tv}(sys::AbstractSystem{Tv}) where Tv  =new(0,0,num_species(sys),sys.grid,sys.physics.data,zeros(Tv,dim_space(sys.grid)))
end

function _fill!(node::BNode{Tv},grid::Grid{Tv},ibnode,ibface) where Tv
    K=grid.bfacenodes[ibnode,ibface]
    node.region=grid.bfaceregions[ibface]
    node.index=K
end





##################################################################
"""
$(TYPEDEF)

Structure holding local node information.

$(TYPEDFIELDS)
"""
mutable struct Node{Tv} <: AbstractGeometryItem 

    """
    Index in grid

    """
    index::Int32
    """
    Inner region number
    """
    region::Int32

    """
    Number of species defined in node
    """
    nspec::Int64

    """
    Number of discretization cell the node is invoked from
    """
    icell::Int64

    """
    Grid
    """
    grid
    
    """
    Physics data
    """
    physics_data


    """
    (deprecated) Coordinates
    """
    coord::Array{Tv,1}
    
    Node{Tv}(sys::AbstractSystem{Tv}) where Tv  =new(0,0,num_species(sys),0, sys.grid,sys.physics.data,zeros(Tv,dim_space(sys.grid)))
end

function _fill!(node::Node{Tv},grid::Grid{Tv},inode,icell) where Tv
    K=cellnode(grid,inode,icell)
    node.region=grid.cellregions[icell]
    node.index=K
    node.icell=icell
end

nodecoord(node::Node)=@views node.grid.coord[:,node.index]

##################################################################
"""
$(TYPEDEF)

Structure holding local edge information.

$(TYPEDFIELDS)
"""
mutable struct Edge{Tv}  <: AbstractGeometryItem

    """
    Index in grid
    """
    index::Int32

    """
    Index 
    """
    node::Vector{Int32}

    """
    Inner region number corresponding to edge
    """
    region::Int32

    """
    Number of species defined in edge
    """
    nspec::Int64

    """
    Number of discretization cell the edge is invoked from
    """
    icell::Int64
    
    """
    Grid
    """
    grid

    """
    Physics data
    """
    physics_data


    """
    (deprecated) Coordinates
    """
    coordK::Array{Tv,1}
    """
    (deprecated) Coordinates
    """
    coordL::Array{Tv,1}

    
    Edge{Tv}(sys::AbstractSystem{Tv}) where Tv  =new(0,[0,0],0,num_species(sys),0,sys.grid,sys.physics.data,zeros(Tv,dim_space(sys.grid)),zeros(Tv,dim_space(sys.grid)))
end


function _fill!(edge::Edge{Tv},grid::Grid{Tv},iedge,icell) where Tv
    K=0
    L=0
    # If we work with projections of fluxes onto edges,
    # we need to ensure that the edges are accessed with the
    # same orientation without regard of the orientation induced
    # by local cell numbering
    if num_edges(grid)>0
        edge.index=celledge(grid,iedge,icell)
        K=grid.edgenodes[1,edge.index]
        L=grid.edgenodes[2,edge.index]
    else
        edge.index=0
        K=celledgenode(grid,1,iedge,icell)
        L=celledgenode(grid,2,iedge,icell)
    end
    edge.region=grid.cellregions[icell]
    edge.node[1]=K
    edge.node[2]=L
    edge.icell=icell
end



nodecoord(edge::Edge,i)=@views edge.grid.coord[:,node.index[i]]

##################################################################
"""
$(TYPEDSIGNATURES)

Return number of species for edge
"""
num_species(edge::Edge{Tv}) where Tv=edge.nspec


##################################################################
"""
$(TYPEDSIGNATURES)
   
Calculate the length of an edge. 
"""
function meas(edge::Edge{Tv}) where Tv
    l::Tv=0.0
    for i=1:dim_space(edge.grid)
        d=edge.grid.coord[i,edge.node[1]]-edge.grid.coord[i,edge.node[2]]
        l=l+d*d
    end
    return sqrt(l)
end




physics_data(item::AbstractGeometryItem)=item.physics_data
####################################################################################

