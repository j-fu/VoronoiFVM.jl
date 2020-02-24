##################################################################
"""
$(TYPEDEF)

Struct holding information for solution array view on subgrid

$(TYPEDFIELDS)
"""
struct SubgridArrayView{Tv,Ta} <: AbstractMatrix{Tv}

    """
    Original array
    """
    sysarray::Ta

    """
    Subgrid for view
    """
    subgrid::SubGrid
end

##################################################################
"""
$(TYPEDSIGNATURES)

Create a view of the solution array on a subgrid.
"""
Base.view(a::AbstractMatrix{Tv},sg::SubGrid) where Tv = SubgridArrayView{Tv,typeof(a)}(a,sg)


##############################################################################
"""
$(TYPEDSIGNATURES)

Accessor method for subgrid array view.
"""
Base.getindex(aview::SubgridArrayView,ispec::Integer,inode::Integer) = aview.sysarray[ispec,aview.subgrid.node_in_parent[inode]]

##############################################################################
"""
$(TYPEDSIGNATURES)

Accessor method for subgrid array view.
"""
@inline function Base.setindex!(aview::SubgridArrayView,v,ispec::Integer,inode::Integer)
    aview.sysarray[ispec,aview.subgrid.node_in_parent[inode]]=v
    return aview
end

##################################################################
"""
$(TYPEDSIGNATURES)
    
Return size of solution array view.
"""
Base.size(a::SubgridArrayView)=(size(a.sysarray,1),size(a.subgrid.node_in_parent,1))


