abstract type AbstractGeometryItem{Tv<:Number, Ti <:Integer} end

##################################################################
"""
$(TYPEDEF)

Structure holding local boundary  node information.

$(TYPEDFIELDS)
"""
mutable struct BNode{Tv, Ti} <: AbstractGeometryItem{Tv, Ti}

    """
    Index in grid
    """
    index::Ti

    """
    Boundary region number
    """
    region::Ti

    """
    Number of species defined in node
    """
    nspec::Ti

    """
    Grid
    """
    grid::Grid{Tv, Ti}
    
    """
    Physics data
    """
    physics_data

    """
    (deprecated) Coordinates
    """
    coord::Array{Tv,1}

    
    BNode{Tv,Ti}(sys::AbstractSystem{Tv,Ti}) where {Tv,Ti}  =new(0,0,num_species(sys),sys.grid,sys.physics.data,zeros(Tv,dim_space(sys.grid)))
end

function _fill!(node::BNode,grid::Grid,ibnode,ibface)
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
mutable struct Node{Tv,Ti} <: AbstractGeometryItem{Tv, Ti} 

    """
    Index in grid

    """
    index::Ti
    """
    Inner region number
    """
    region::Ti

    """
    Number of species defined in node
    """
    nspec::Ti

    """
    Number of discretization cell the node is invoked from
    """
    icell::Ti

    """
    Grid
    """
    grid::Grid{Tv,Ti}
    
    """
    Physics data
    """
    physics_data


    """
    (deprecated) Coordinates
    """
    coord::Array{Tv,1}
    
    Node{Tv,Ti}(sys::AbstractSystem{Tv,Ti}) where {Tv,Ti}  =new(zero(Ti),0,num_species(sys),0, sys.grid,sys.physics.data,zeros(Tv,dim_space(sys.grid)))
end

function _fill!(node::Node,grid::Grid,inode,icell)
    K=cellnode(grid,inode,icell)
    node.region=grid.cellregions[icell]
    node.index=K
    node.icell=icell
end


Base.size(node::Node)=(size(node.grid.coord)[1],)

Base.getindex(node::Node, idim)= node.grid.coord[idim,node.index]

##################################################################
"""
$(TYPEDEF)

Structure holding local edge information.

$(TYPEDFIELDS)
"""
mutable struct Edge{Tv,Ti}  <: AbstractGeometryItem{Tv, Ti}

    """
    Index in grid
    """
    index::Ti

    """
    Index 
    """
    node::Vector{Ti}

    """
    Inner region number corresponding to edge
    """
    region::Ti

    """
    Number of species defined in edge
    """
    nspec::Ti

    """
    Number of discretization cell the edge is invoked from
    """
    icell::Ti
    
    """
    Grid
    """
    grid::Grid{Tv, Ti}

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

    
    Edge{Tv,Ti}(sys::AbstractSystem{Tv,Ti}) where {Tv,Ti}  =new(0,[0,0],0,num_species(sys),0,sys.grid,sys.physics.data,zeros(Tv,dim_space(sys.grid)),zeros(Tv,dim_space(sys.grid)))
end


function _fill!(edge::Edge,grid::Grid,iedge,icell) where Tv
    if num_edges(grid)>0
        # If we work with projections of fluxes onto edges,
        # we need to ensure that the edges are accessed with the
        # same orientation without regard of the orientation induced
        # by local cell numbering
        edge.index=celledge(grid,iedge,icell)
        edge.node[1]=grid.edgenodes[1,edge.index]
        edge.node[2]=grid.edgenodes[2,edge.index]
    else
        edge.index=0
        edge.node[1]=celledgenode(grid,1,iedge,icell)
        edge.node[2]=celledgenode(grid,2,iedge,icell)
    end
    edge.region=grid.cellregions[icell]
    edge.icell=icell
end


Base.size(edge::Edge)=(size(edge.grid.coord)[1],2)
Base.getindex(edge::Edge, idim,inode)= edge.grid.coord[idim,edge.node[inode]]

##################################################################
"""
$(TYPEDSIGNATURES)

Return number of species for edge
"""
num_species(edge::Edge)=edge.nspec


##################################################################
"""
$(TYPEDSIGNATURES)
   
Calculate the length of an edge. 
"""
function meas(edge::Edge)
    l=0.0
    for i=1:dim_space(edge.grid)
        d=edge.grid.coord[i,edge.node[1]]-edge.grid.coord[i,edge.node[2]]
        l=l+d*d
    end
    return sqrt(l)
end




physics_data(item::AbstractGeometryItem)=item.physics_data
####################################################################################

