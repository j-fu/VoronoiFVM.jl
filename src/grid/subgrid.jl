#=
Definition of subgrid data type.
=#

##################################################################
"""
$(TYPEDEF)
    
Subgrid of parent grid (mainly for visualization purposes). Intended
to hold support of species which are not defined everywhere.

$(TYPEDFIELDS)
"""
struct SubGrid{Tc} <: AbstractGrid


    """
    Parent Grid
    """
    parent::Grid

    
    """
    Incidence between subgrid node numbers and node numbers
    in parent.
    """
    node_in_parent::Array{Int32,1}

    
    """ 
    2D Array of coordinates per grid node
    """
    coord::Array{Tc,2}

    
    """
    2D Array of node numbers per grid cell
    """
    cellnodes::Array{Int32,2}

    
end

function Base.show(io::IO,subgrid::SubGrid)
    str=@sprintf("%s(dim_space=%d, num_nodes=%d, num_cells=%d)",
                 typeof(subgrid),dim_space(subgrid),num_nodes(subgrid), num_cells(subgrid))
    println(io,str)
end

##################################################################
# Default transform for subgrid creation
function _copytransform!(a::AbstractArray,b::AbstractArray)
    for i=1:length(a)
        a[i]=b[i]
    end
end

##################################################################
"""
$(TYPEDSIGNATURES)

Create subgrid of list of regions.
"""
function subgrid(parent::Grid,
                 subregions::AbstractArray;
                 transform::Function=_copytransform!,
                 boundary=false)
    Tc=Base.eltype(parent)
    
    @inline function insubregions(xreg)
        for i in eachindex(subregions)
            if subregions[i]==xreg
                return true
            end
        end
        return false
    end

    
    if boundary
        xregions=parent.bfaceregions
        xnodes=parent.bfacenodes
        sub_gdim=dim_grid(parent)-1
    else
        xregions=parent.cellregions
        xnodes=parent.cellnodes
        sub_gdim=dim_grid(parent)
    end
    
    nodemark=zeros(Int32,num_nodes(parent))
    ncn=size(xnodes,1)
    
    nsubcells=0
    nsubnodes=0
    for icell in eachindex(xregions)
        if insubregions(xregions[icell])
            nsubcells+=1
            for inode=1:ncn
                ipnode=xnodes[inode,icell]
                if nodemark[ipnode]==0
                    nsubnodes+=1
                    nodemark[ipnode]=nsubnodes
                end
            end
        end
    end
    
    sub_cellnodes=zeros(Int32,ncn,nsubcells)
    sub_nip=zeros(Int32,nsubnodes)
    for inode in eachindex(nodemark)
        if nodemark[inode]>0
            sub_nip[nodemark[inode]]=inode
        end
    end
    
    isubcell=0
    for icell in eachindex(xregions)
        if insubregions(xregions[icell])
            isubcell+=1
            for inode=1:ncn
                ipnode=xnodes[inode,icell]
                sub_cellnodes[inode,isubcell]=nodemark[ipnode]
            end
        end
    end

    localcoord=zeros(Tc,sub_gdim,nsubnodes)
    @views for inode=1:nsubnodes
        transform(localcoord[:,inode],parent.coord[:,sub_nip[inode]])
    end
    
    return SubGrid(parent,sub_nip,localcoord,sub_cellnodes)
end


