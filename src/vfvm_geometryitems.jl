"""
   $(TYPEDEF)

Abstract type for geometry items (node,bnode,edge, bedge)
"""
abstract type AbstractGeometryItem{Tv<:Number, Ti <:Integer} end


time(item::AbstractGeometryItem)=item.time
embedparam(item::AbstractGeometryItem)=item.embedparam
region(item::AbstractGeometryItem)=item.region


"""
   $(TYPEDEF)

Abstract type for nodes. 

`node[idim]` gives the the corresponding coordinate.
"""
abstract type AbstractNode{Tv<:Number, Ti <:Integer}  <: AbstractGeometryItem{Tv, Ti} end
Base.size(node::AbstractNode)=(size(node.coord)[1],)
Base.getindex(node::AbstractNode, idim)=@inbounds  node.coord[idim,node.index]


"""
    $(TYPEDEF)

Abstract type for data on nodes.
`u[ispec]` accesses value of species at this node.
"""
abstract type AbstractNodeData{T<: Number} <: AbstractVector{T} end
Base.size(u::AbstractNodeData)=(u.nspec,1)
Base.getindex(u::AbstractNodeData,i)=@inbounds u.val[i]
Base.setindex!(f::AbstractNodeData,v,i)=@inbounds f.val[i]=v


struct DParameters{T<:Number} <: AbstractVector{T}
    val::Vector{T}
    offset::Int32
end

Base.size(p::DParameters)=(length(p.val)-p.offset,1)
Base.getindex(p::DParameters,i)=@inbounds p.val[p.offset+i]

function parameters(u::AbstractNodeData{T}) where {T<:Number}
    DParameters(u.val,u.nspec)
end



"""
   $(TYPEDEF)

Abstract type for edges 

`edge[idim,inode]` gives coordinate of node.
"""
abstract type AbstractEdge{Tv<:Number, Ti <:Integer}  <: AbstractGeometryItem{Tv, Ti} end
Base.size(edge::AbstractEdge)=(size(edge.coord)[1],2)
Base.getindex(edge::AbstractEdge, idim,inode)=@inbounds  edge.coord[idim,edge.node[inode]]

"""
    $(TYPEDEF)

Abstract type for data on edges.
`u[ispec,inode]` accesses value of species at corresponding node.
"""
abstract type AbstractEdgeData{T<: Number} <: AbstractMatrix{T} end
Base.size(u::AbstractEdgeData)=(u.n1,2)
Base.getindex(u::AbstractEdgeData,i,j)=@inbounds u.val[(j-1)*u.n1+i]

function parameters(u::AbstractEdgeData{T}) where {T<:Number}
    DParameters(u.val,2*u.n1)
end


##################################################################
"""
$(TYPEDEF)

Structure holding local node information.

$(TYPEDFIELDS)
"""
mutable struct Node{Tv,Ti} <: AbstractNode{Tv, Ti} 

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


    """
    Grid cell nodes
    """
    cellnodes::Array{Ti,2}

    
    """
    Grid cell regions
    """
    cellregions::Vector{Ti}

    """
    System time
    """
    time::Tv
    
    """
    Current value of embedding parameter
    """
    embedparam::Tv
    
    Node{Tv,Ti}(sys::AbstractSystem{Tv,Ti}) where {Tv,Ti}  =new(zero(Ti),0,
                                                                num_species(sys),0,
                                                                coordinates(sys.grid),
                                                                sys.grid[CellNodes],
                                                                sys.grid[CellRegions],
                                                                0,0
                                                                )
end


@inline function _fill!(node::Node,inode,icell)
    node.region=node.cellregions[icell]
    node.index=node.cellnodes[inode,icell]
    node.icell=icell
    nothing
end


"""
    $(TYPEDEF)

Unknown data on node. 
"""
struct NodeUnknowns{T,Tv,Ti} <:AbstractNodeData{T} 
    val::Vector{T}
    nspec::Ti
    geom::Node{Tv,Ti}
end

@inline unknowns(node::Node{Tv,Ti},u::Vector{T}) where {T,Tv,Ti} = NodeUnknowns{T,Tv,Ti}(u,node.nspec,node)

"""
    $(TYPEDEF)

RHS data on node. 
"""
struct NodeRHS{T,Tv,Ti} <:AbstractNodeData{T}
    val::Vector{T}
    nspec::Ti
    geom::Node{Tv,Ti}
end

@inline rhs(node::Node{Tv,Ti}, f::Vector{T}) where {T,Tv,Ti} = NodeRHS{T,Tv,Ti}(f,node.nspec,node)


##################################################################
"""
$(TYPEDEF)

Structure holding local boundary  node information.

$(TYPEDFIELDS)
"""
mutable struct BNode{Tv, Ti} <: AbstractNode{Tv, Ti}

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
    
    Dirichlet::Tv

    """
    System time
    """
    time::Tv

    """
    Current value of embedding parameter
    """
    embedparam::Tv

    dirichlet_value::Vector{Tv}
    
    BNode{Tv,Ti}(sys::AbstractSystem{Tv,Ti}) where {Tv,Ti}  =new(0,0,0,0,zeros(Ti,2),
                                                                 num_species(sys),
                                                                 coordinates(sys.grid),
                                                                 sys.grid[BFaceNodes],
                                                                 sys.grid[BFaceRegions],
                                                                 sys.grid[CellRegions],
                                                                 sys.grid[BFaceCells],
                                                                 Dirichlet,0,0,
                                                                 zeros(Tv,num_species(sys))
                                                                 )
end


@inline function _fill0!(node::BNode,ibnode,ibface)
    node.ibface=ibface
    node.ibnode=ibnode
    node.region=node.bfaceregions[ibface]
    node.index=node.bfacenodes[ibnode,ibface]
    nothing
end


@inline function _fill!(node::BNode,ibnode,ibface)
    _fill0!(node,ibnode,ibface)
    node.cellregions[1]=0
    node.cellregions[2]=0
    for i=1:num_targets(node.bfacecells,ibface)
        icell=node.bfacecells[i,ibface]
        node.cellregions[i]=node.allcellregions[icell]
    end
end




struct BNodeUnknowns{T,Tv,Ti} <:AbstractNodeData{T} 
    val::Vector{T}
    nspec::Ti
    geom::BNode{Tv,Ti}
end

@inline unknowns(bnode::BNode{Tv,Ti},u::Vector{T}) where {T,Tv,Ti} = BNodeUnknowns{T,Tv,Ti}(u,bnode.nspec,bnode)


struct BNodeRHS{T,Tv,Ti} <:AbstractNodeData{T} 
    val::Vector{T}
    nspec::Ti
    geom::BNode{Tv,Ti}
end

@inline rhs(bnode::BNode{Tv,Ti}, f::Vector{T}) where {T,Tv,Ti} = BNodeRHS{T,Tv,Ti}(f,bnode.nspec,bnode)




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

    """
    System time
    """
    time::Tv

    """
    Current value of embedding parameter
    """
    embedparam::Tv
    
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
    edge.time=0
    edge.embedparam=0
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


struct EdgeUnknowns{T,Tv,Ti} <:AbstractEdgeData{T} 
    val::Vector{T}
    n1::Ti
    geom::Edge{Tv,Ti}
end

@inline unknowns(edge::Edge{Tv,Ti},u::Vector{T}) where {T,Tv,Ti} = EdgeUnknowns{T,Tv,Ti}(u,edge.nspec,edge)


struct EdgeRHS{T,Tv,Ti} <:AbstractNodeData{T} 
    val::Vector{T}
    nspec::Ti
    geom::Edge{Tv,Ti}
end

@inline rhs(edge::Edge{Tv,Ti}, f::Vector{T}) where {T,Tv,Ti} = EdgeRHS{T,Tv,Ti}(f,edge.nspec,edge)








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

    """
    System time
    """
    time::Tv

    """
    Current value of embedding parameter
    """
    embedparam::Tv

    
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
    bedge.time=0
    bedge.embedparam=0
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

struct BEdgeUnknowns{T,Tv,Ti} <:AbstractEdgeData{T} 
    val::Vector{T}
    n1::Ti
    geom::BEdge{Tv,Ti}
end

@inline unknowns(edge::BEdge{Tv,Ti},u::Vector{T}) where {T,Tv,Ti} = BEdgeUnknowns{T,Tv,Ti}(u,edge.nspec,edge)


struct BEdgeRHS{T,Tv,Ti} <:AbstractNodeData{T} 
    val::Vector{T}
    nspec::Ti
    geom::BEdge{Tv,Ti}
end

@inline rhs(edge::BEdge{Tv,Ti}, f::Vector{T}) where {T,Tv,Ti}= BEdgeRHS{T,Tv,Ti}(f,edge.nspec,edge)


##################################################################
"""
$(TYPEDSIGNATURES)

Return number of species for edge
"""
@inline num_species(edge::AbstractEdge)=edge.nspec


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






###############################################################
# Deprecation warnings here ?
"""
$(TYPEDEF)

Wrapper struct for viewing unknowns passed to callback functions
    
$(TYPEDFIELDS)
"""
struct VectorUnknowns{T} <:AbstractVector{T} 
    val::Vector{T}
    n::Int64
    offset::Int64
end


"""
$(TYPEDSIGNATURES)

Construct vector unknowns on edge.
"""
unknowns(edge::AbstractEdge, u::Vector{T},i) where T = VectorUnknowns{T}(u,edge.nspec,(i-1)*(edge.nspec))
Base.getindex(u::VectorUnknowns,i)=@inbounds u.val[u.offset+i]
Base.size(u::VectorUnknowns)=(u.n,)



# For backward compatibility
unknowns(edge,u::AbstractEdgeData)=u
# For backward compatibility
unknowns(edge::Edge, u::AbstractEdgeData{T},i) where T = VectorUnknowns{T}(u.val,edge.nspec,(i-1)*(edge.nspec))

"""
$(TYPEDSIGNATURES)

Solution view on first edge node
"""
viewK(edge::AbstractEdge,u)=unknowns(edge,u,1)


"""
$(TYPEDSIGNATURES)

Solution view on second edge node
"""
viewL(edge::AbstractEdge,u)=unknowns(edge,u,2)

