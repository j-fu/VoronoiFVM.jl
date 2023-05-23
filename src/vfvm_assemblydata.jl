abstract type AbstractAssemblyData{Tv,Ti} end


struct CellWiseAssemblyData{Tv,Ti} <: AbstractAssemblyData{Tv,Ti}
    """
    Precomputed geometry factors for cell nodes
    """
    cellnodefactors::Array{Tv, 2}

    """
    Precomputed geometry factors for cell edges
    """
    celledgefactors::Array{Tv, 2}

    """
    Cell regions
    """
    cellregions::Array{Ti,1}
end


struct EdgeWiseAssemblyData{Tv,Ti}<:AbstractAssemblyData{Tv,Ti}
    """
        Precomputed geometry factors for cell nodes
    """
    nodefactors::SparseMatrixCSC{Tv, Ti}

    """
    Precomputed geometry factors for cell edges
    """
    edgefactors::SparseMatrixCSC{Tv, Ti}
end



noderange(asmdata::CellWiseAssemblyData)=1:size(asmdata.cellnodefactors,2)
edgerange(asmdata::CellWiseAssemblyData)=1:size(asmdata.celledgefactors,2)

noderange(asmdata::CellWiseAssemblyData,icell)=1:size(asmdata.cellnodefactors,1)
edgerange(asmdata::CellWiseAssemblyData,icell)=1:size(asmdata.celledgefactors,1)



function _fill!(node::Node,asmdata::CellWiseAssemblyData,inode,icell)
    node.fac=asmdata.cellnodefactors[inode,icell]
    _fill!(node, inode, icell)
end

function _fill!(edge::Edge,asmdata::CellWiseAssemblyData,iedge,icell)
    edge.fac=asmdata.celledgefactors[iedge,icell]
    _fill!(edge, iedge, icell)
end


noderange(asmdata::EdgeWiseAssemblyData)=1:size(asmdata.nodefactors,2)
edgerange(asmdata::EdgeWiseAssemblyData)=1:size(asmdata.edgefactors,2)

noderange(asmdata::EdgeWiseAssemblyData, inode)=nzrange(asmdata.nodefactors, inode)
edgerange(asmdata::EdgeWiseAssemblyData, iedge)=nzrange(asmdata.edgefactors, iedge)


function _fill!(node::Node,asmdata::EdgeWiseAssemblyData,k,inode)
    _xfill!(node, asmdata.nodefactors.rowval[k], inode)
    node.fac = asmdata.nodefactors.nzval[k]
end

function _fill!(edge::Edge,asmdata::EdgeWiseAssemblyData,k,iedge)
    _xfill!(edge, asmdata.edgefactors.rowval[k], iedge)
    edge.fac = asmdata.edgefactors.nzval[k]
end


