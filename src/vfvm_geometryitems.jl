##################################################################
"""
$(TYPEDEF)

Structure holding local boundary  node information.

$(TYPEDFIELDS)
"""
mutable struct BNode{Tv}

    """
    Index in grid
    """
    index::Int32

    """
    Boundary region number
    """
    region::Int32

    """
    1D Array of node coordinates
    """
    coord::Array{Tv,1}

    """
    Number of species defined in node
    """
    nspec::Int64
    BNode{Tv}(sys::AbstractSystem{Tv}) where Tv  =new(0,0,zeros(Tv,dim_space(sys.grid)),num_species(sys))
end

function _fill!(node::BNode{Tv},grid::Grid{Tv},ibnode,ibface) where Tv
    K=grid.bfacenodes[ibnode,ibface]
    node.region=grid.bfaceregions[ibface]
    node.index=K
    for i=1:length(node.coord)
        node.coord[i]=grid.coord[i,K]
    end
end





##################################################################
"""
$(TYPEDEF)

Structure holding local node information.

$(TYPEDFIELDS)
"""
mutable struct Node{Tv}

    """
    Index in grid

    """
    index::Int32
    """
    Inner region number
    """
    region::Int32

    """
    1D Array of node coordinates
    """
    coord::Array{Tv,1}

    """
    Number of species defined in node
    """
    nspec::Int64

    """
    Number of discretization cell the node is invoked from
    """
    icell::Int64

    Node{Tv}(sys::AbstractSystem{Tv}) where Tv  =new(0,0,zeros(Tv,dim_space(sys.grid)),num_species(sys),0)
end

function _fill!(node::Node{Tv},grid::Grid{Tv},inode,icell) where Tv
    K=cellnode(grid,inode,icell)
    node.region=grid.cellregions[icell]
    node.index=K
    node.icell=icell
    for i=1:length(node.coord)
        node.coord[i]=grid.coord[i,K]
    end
end

##################################################################
"""
$(TYPEDEF)

Structure holding local edge information.

$(TYPEDFIELDS)
"""
mutable struct Edge{Tv}

    """
    Index in grid
    """
    index::Int32

    """
    Index of first node
    """
    nodeK::Int32

    """
    Index of second node
    """
    nodeL::Int32

    """
    Inner region number corresponding to edge
    """
    region::Int32

    """
    1D Array of first node coordinates
    """
    coordK::Array{Tv,1}

    """
    1D Array of second node coordinates
    """
    coordL::Array{Tv,1}

    """
    Number of species defined in edge
    """
    nspec::Int64

    """
    Number of discretization cell the edge is invoked from
    """
    icell::Int64

    Edge{Tv}(sys::AbstractSystem{Tv}) where Tv  =new(0,0,0,0,zeros(Tv,dim_space(sys.grid)),zeros(Tv,dim_space(sys.grid)),num_species(sys),0)
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
    edge.nodeK=K
    edge.nodeL=L
    edge.icell=icell
    @. @inbounds @views edge.coordK=grid.coord[:,K]
    @. @inbounds @views edge.coordL=grid.coord[:,L]
end


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
function edgelength(edge::Edge{Tv}) where Tv
    l::Tv=0.0
    for i=1:length(edge.coordK)
        d=edge.coordK[i]-edge.coordL[i]
        l=l+d*d
    end
    return sqrt(l)
end

"""
$(TYPEDSIGNATURES)

Solution view on first edge node
"""
@inline viewK(edge::Edge{Tv},u::AbstractArray) where Tv=@views u[1:edge.nspec]


"""
$(TYPEDSIGNATURES)

Solution view on second edge node
"""
@inline viewL(edge::Edge{Tv},u::AbstractArray) where Tv=@views u[edge.nspec+1:2*edge.nspec]

"""
$(TYPEDSIGNATURES)

Solution view on first edge node
"""
@inline viewK(nspec::Int64,u::AbstractArray)=@views u[1:nspec]


"""
$(TYPEDSIGNATURES)

Solution view on second edge node
"""
@inline viewL(nspec::Int64,u::AbstractArray)=@views u[nspec+1:2*nspec]



function Base.show(io::IO,sys::AbstractSystem) where Tc
    str=@sprintf("%s(num_species=%d)",typeof(sys),sys.physics.num_species)
    println(io,str)
end


####################################################################################

