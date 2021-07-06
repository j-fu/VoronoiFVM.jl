abstract type AbstractGeometryItem{Tv<:Number, Ti <:Integer} end

abstract type AbstractEdge{Tv<:Number, Ti <:Integer}  <: AbstractGeometryItem{Tv, Ti} end

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
    BFace number it is called from
    """
    ibface::Ti

    
    """
    local node number
    """
    ibnode::Ti

    """
    Boundary region number
    """
    region::Ti


    cellregions::Vector{Ti}

    """
    Number of species defined in node
    """
    nspec::Ti

    """
    Grid coordinates
    """
    coord::Matrix{Tv}


    bfacenodes::Array{Ti,2}

    bfaceregions::Vector{Ti}

    allcellregions::Vector{Ti}
    
    bfacecells::ExtendableGrids.Adjacency{Ti}
    
    
    BNode{Tv,Ti}(sys::AbstractSystem{Tv,Ti}) where {Tv,Ti}  =new(0,0,0,0,zeros(Ti,2),
                                                                 num_species(sys),
                                                                 coordinates(sys.grid),
                                                                 sys.grid[BFaceNodes],
                                                                 sys.grid[BFaceRegions],
                                                                 sys.grid[CellRegions],
                                                                 sys.grid[BFaceCells],
                                                                 )
end


@inline function _fill0!(node::BNode,ibnode,ibface)
    node.ibface=ibface
    node.ibnode=ibnode
    node.region=node.bfaceregions[ibface]
    node.index=node.bfacenodes[ibnode,ibface]
    nothing
end


function _fill!(node::BNode,ibnode,ibface)
    _fill0!(node,ibnode,ibface)
    node.cellregions[1]=0
    node.cellregions[2]=0
    for i=1:num_targets(node.bfacecells,ibface)
        icell=node.bfacecells[i,ibface]
        node.cellregions[i]=node.allcellregions[icell]
    end
end





Base.size(bnode::BNode)=(size(bnode.coord)[1],)
Base.getindex(bnode::BNode, idim)= bnode.coord[idim,bnode.index]




struct BNodeUnknowns{T} <:AbstractMatrix{T} 
    u::Vector{T}
    bnode::BNode
end

unknowns(bnode::BNode,u::Vector{T}) where T = BNodeUnknowns{T}(u,bnode)
Base.getindex(u::BNodeUnknowns,i)=@inbounds u.u[i]


struct BNodeRHS{T} <:AbstractVector{T}
    f::Vector{T}
    bnode::BNode
end

Base.size(f::BNodeRHS)=size(f.f)


rhs(bnode::BNode, f::Vector{T}) where T = BNodeRHS(f,bnode)
Base.getindex(f::BNodeRHS,i)=@inbounds f.f[i]




Base.setindex!(f::BNodeRHS,v,i)=@inbounds f.f[i]=v



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
    Grid coordinates
    """
    coord::Matrix{Tv}


    cellnodes::Array{Ti,2}

    cellregions::Vector{Ti}

    
    Node{Tv,Ti}(sys::AbstractSystem{Tv,Ti}) where {Tv,Ti}  =new(zero(Ti),0,
                                                                num_species(sys),0,
                                                                coordinates(sys.grid),
                                                                sys.grid[CellNodes],
                                                                sys.grid[CellRegions]
                                                                )
end


@inline function _fill!(node::Node,inode,icell)
    node.region=node.cellregions[icell]
    node.index=node.cellnodes[inode,icell]
    node.icell=icell
    nothing
end


Base.size(node::Node)=(size(node.coord)[1],)

Base.getindex(node::Node, idim)= node.coord[idim,node.index]

##################################################################
"""
$(TYPEDEF)

Structure holding local edge information.

$(TYPEDFIELDS)
"""
mutable struct Edge{Tv,Ti}  <: AbstractEdge{Tv, Ti}

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
    Grid coordinates
    """
    coord::Matrix{Tv}

    
    cellx::Array{Ti,2}
    edgenodes::Array{Ti,2}
    cellregions::Vector{Ti}
    has_celledges::Bool
    
    Edge{Tv,Ti}(::Nothing) where {Tv,Ti}  =new()
end


function  Edge{Tv,Ti}(sys::AbstractSystem{Tv,Ti}) where {Tv,Ti} 
    edge=Edge{Tv,Ti}(nothing)

    edge.index=0
    edge.node=[0,0]
    edge.region=0
    edge.nspec=num_species(sys)
    edge.icell=0
    edge.coord=coordinates(sys.grid)
    geom=sys.grid[CellGeometries][1]
    if haskey(sys.grid,CellEdges)
        edge.cellx=sys.grid[CellEdges]
        edge.edgenodes=sys.grid[EdgeNodes]
        edge.has_celledges=true
    else
        edge.cellx=sys.grid[CellNodes]
        edge.edgenodes=local_celledgenodes(geom)
        edge.has_celledges=false
    end
    edge.cellregions=sys.grid[CellRegions]
    edge
end


@inline function _fill!(edge::Edge,iedge,icell)
    if edge.has_celledges #  cellx==celledges, edgenodes==global_edgenodes
        # If we work with projections of fluxes onto edges,
        # we need to ensure that the edges are accessed with the
        # same orientation without regard of the orientation induced
        # by local cell numbering
        edge.index=edge.cellx[iedge,icell]
        edge.node[1]=edge.edgenodes[1,edge.index]
        edge.node[2]=edge.edgenodes[2,edge.index]
    else # cx==cellnodes, edgenodes== local_edgenodes
        edge.index=0
        edge.node[1]=edge.cellx[edge.edgenodes[1,iedge],icell]
        edge.node[2]=edge.cellx[edge.edgenodes[2,iedge],icell]
    end
    edge.region=edge.cellregions[icell]
    edge.icell=icell
    nothing
end

function _fill!(edge::Edge,cellx,edgenodes,cellregions,iedge,icell, has_celledges) where Tv
    if has_celledges #  cellx==celledges, edgenodes==global_edgenodes
        # If we work with projections of fluxes onto edges,
        # we need to ensure that the edges are accessed with the
        # same orientation without regard of the orientation induced
        # by local cell numbering
        edge.index=cellx[iedge,icell]
        edge.node[1]=edgenodes[1,edge.index]
        edge.node[2]=edgenodes[2,edge.index]
    else # cx==cellnodes, edgenodes== local_edgenodes
        edge.index=0
        edge.node[1]=cellx[edgenodes[1,iedge],icell]
        edge.node[2]=cellx[edgenodes[2,iedge],icell]
    end
    edge.region=cellregions[icell]
    edge.icell=icell
end


Base.size(edge::Edge)=(size(edge.coord)[1],2)
Base.getindex(edge::Edge, idim,inode)= edge.coord[idim,edge.node[inode]]

struct EdgeUnknowns{T} <:AbstractMatrix{T} 
    u::Vector{T}
    n1::Int64
    edge::AbstractEdge
end

Base.size(u::EdgeUnknowns)=(u.n1,2)

unknowns(edge::AbstractEdge,u::Vector{T}) where T = EdgeUnknowns{T}(u,edge.nspec,edge)
Base.getindex(u::EdgeUnknowns,i,j)=@inbounds u.u[(j-1)*u.n1+i]


struct EdgeRHS{T} <:AbstractVector{T} 
    f::Vector{T}
    edge::Edge
end


rhs(edge::Edge, f::Vector{T}) where T = EdgeRHS(f,edge)
Base.getindex(f::EdgeRHS,i)=@inbounds f.f[i]

Base.setindex!(f::EdgeRHS,v,i)=@inbounds f.f[i]=v








##################################################################
"""
$(TYPEDEF)

Structure holding local edge information.

$(TYPEDFIELDS)
"""
mutable struct BEdge{Tv,Ti}  <: AbstractEdge{Tv, Ti}

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
    Grid coordinates
    """
    coord::Matrix{Tv}

    bedgenodes::Array{Ti,2}
    bfaceedges::Array{Ti,2}
    bfaceregions::Vector{Ti}
    
    BEdge{Tv,Ti}(::Nothing) where {Tv,Ti}  =new()
end


function  BEdge{Tv,Ti}(sys::AbstractSystem{Tv,Ti}) where {Tv,Ti} 
    bedge=BEdge{Tv,Ti}(nothing)

    bedge.index=0
    bedge.node=[0,0]
    bedge.region=0
    bedge.nspec=num_species(sys)
    bedge.icell=0
    bedge.coord=coordinates(sys.grid)
       
    bedge.bedgenodes=sys.grid[BEdgeNodes]
    bedge.bfaceedges=sys.grid[BFaceEdges]
    bedge.bfaceregions=sys.grid[BFaceRegions]
    bedge
end


@inline function _fill!(bedge::BEdge,ibedge,ibface)

    bedge.index   = bedge.bfaceedges[ibedge, ibface]
    bedge.node[1] = bedge.bedgenodes[1, bedge.index]
    bedge.node[2] = bedge.bedgenodes[2, bedge.index]
    
    bedge.region = bedge.bfaceregions[ibface]
    bedge.icell  = ibface

    nothing
end

Base.size(bedge::BEdge)=(size(edge.coord)[1],2)
Base.getindex(bedge::BEdge, idim,inode)= bedge.coord[idim,bedge.node[inode]]

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
function meas(edge::AbstractEdge)
    l=0.0
    for i=1:size(edge.coord)[1]
        d=edge.coord[i,edge.node[1]]-edge.coord[i,edge.node[2]]
        l=l+d*d
    end
    return sqrt(l)
end

edgelength(edge::AbstractEdge)=meas(edge)


function project(edge::Edge,vec)
    vh=0.0
    for i=1:size(edge.coord)[1]
        vh+=(edge.coord[i,edge.node[2]]-edge.coord[i,edge.node[1]])*vec[i]
    end
    return vh
end






##################################################################
"""
$(TYPEDEF)

Wrapper struct for viewing unknowns passed to callback functions
    
$(TYPEDFIELDS)
"""
struct VectorUnknowns{T} <:AbstractVector{T} 
    u::Vector{T}
    n::Int64
    offset::Int64
end


"""
$(TYPEDSIGNATURES)

Construct vector unknowns on edge.
"""
unknowns(edge::Edge, u::Vector{T},i) where T = VectorUnknowns{T}(u,edge.nspec,(i-1)*(edge.nspec))
Base.getindex(u::VectorUnknowns,i)=@inbounds u.u[u.offset+i]
Base.size(u::VectorUnknowns)=(u.n,)


"""
$(TYPEDSIGNATURES)

Construct vector unknowns on node
"""
unknowns(node::Node, u)=u

"""
$(TYPEDSIGNATURES)

Construct vector unknowns on bnode
"""
unknowns(node::BNode, u)=u


################################################
# Deprecation warnings here!

# For backward compatibility
unknowns(edge,u::EdgeUnknowns)=u
# For backward compatibility
unknowns(edge::Edge, u::EdgeUnknowns{T},i) where T = VectorUnknowns{T}(u.u,edge.nspec,(i-1)*(edge.nspec))

"""
$(TYPEDSIGNATURES)

Solution view on first edge node
"""
viewK(edge::Edge,u)=unknowns(edge,u,1)


"""
$(TYPEDSIGNATURES)

Solution view on second edge node
"""
viewL(edge::Edge,u)=unknowns(edge,u,2)

