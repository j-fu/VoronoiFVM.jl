"""
    $(TYPEDEF)

Assembly of residual and Jacobian comes in two flavors, cellwise and edgewise assembly loops, see [`VoronoiFVM.System(grid;kwargs...)`](@ref).
The necessary data for assembly are held in structs which are subtypes of `AbstractAssemblyData`.
"""
abstract type AbstractAssemblyData{Tv, Ti} end

"""
    $(TYPEDEF)
Data for cellwise assembly.


$(TYPEDFIELDS)
"""
struct CellwiseAssemblyData{Tv, Ti} <: AbstractAssemblyData{Tv, Ti}
    """
        Precomputed geometry factors for cell nodes.
        This is a `ncells x nnodes_per_cell` full matrix.
    """
    nodefactors::Array{Tv,2}

    """
        Precomputed geometry factors for cell edges
        This is a `ncells x nedge_per_cell` full matrix.
    """
    edgefactors::Array{Tv,2}

end


"""
    $(TYPEDEF)


$(TYPEDFIELDS)
"""
struct EdgewiseAssemblyData{Tv,Ti} <: AbstractAssemblyData{Tv,Ti}
    """
        Precomputed geometry factors for  nodes.
        This is a `nnodes x nregions` sparse matrix.
    """
    nodefactors::SparseMatrixCSC{Tv,Ti}

    """
    Precomputed geometry factors for  edges
        This is a `nedges x nregions` sparse matrix.
    """
    edgefactors::SparseMatrixCSC{Tv,Ti}
end

"""
    nodebatch(assemblydata)

Outer range for node assembly loop. 
"""
function nodebatch end

"""
    noderange(assemblydata, i)

Inner range for node assembly loop. 
"""
function noderange end

"""
    nodebatch(assemblydata)

Outer range for edge assembly loop. 
"""
function edgebatch end

"""
    edgerange(assemblydata, i)

Inner range for edge assembly loop. 
"""
function edgerange end




nodebatch(asmdata::CellwiseAssemblyData{Tv,Ti}) where {Tv,Ti} =
    1:size(asmdata.nodefactors, 2)
edgebatch(asmdata::CellwiseAssemblyData{Tv,Ti}) where {Tv,Ti} =
    1:size(asmdata.edgefactors, 2)

nodebatch(asmdata::EdgewiseAssemblyData{Tv,Ti}) where {Tv,Ti} =
    1:size(asmdata.nodefactors, 2)
edgebatch(asmdata::EdgewiseAssemblyData{Tv,Ti}) where {Tv,Ti} =
    1:size(asmdata.edgefactors, 2)


noderange(asmdata::EdgewiseAssemblyData{Tv,Ti}, inode) where {Tv,Ti} =
    nzrange(asmdata.nodefactors, inode)
edgerange(asmdata::EdgewiseAssemblyData{Tv,Ti}, iedge) where {Tv,Ti} =
    nzrange(asmdata.edgefactors, iedge)

noderange(asmdata::CellwiseAssemblyData{Tv,Ti}, inode) where {Tv,Ti} =
    1:size(asmdata.nodefactors, 1)
edgerange(asmdata::CellwiseAssemblyData{Tv,Ti}, iedge) where {Tv,Ti} =
    1:size(asmdata.edgefactors, 1)

"""
    $(SIGNATURES)

Fill node with the help of assemblydata.
"""
function _fill!(
    node::Node,
    asmdata::CellwiseAssemblyData{Tv,Ti},
    inode,
    icell,
) where {Tv,Ti}
    node.index = node.cellnodes[inode, icell]
    node.region = node.cellregions[icell]
    node.fac = asmdata.nodefactors[inode, icell]
    node.icell = icell
end

"""
    $(SIGNATURES)

Fill boundary node with the help of assemblydata.
"""
function _fill!(
    node::BNode,
    asmdata::CellwiseAssemblyData{Tv,Ti},
    ibnode,
    ibface,
) where {Tv,Ti}
    node.ibface = ibface
    node.ibnode = ibnode
    node.region = node.bfaceregions[ibface]
    node.index = node.bfacenodes[ibnode, ibface]
    node.cellregions[1] = 0
    node.cellregions[2] = 0
    for i = 1:num_targets(node.bfacecells, ibface)
        icell = node.bfacecells[i, ibface]
        node.cellregions[i] = node.allcellregions[icell]
    end
    node.fac = asmdata.nodefactors[ibnode, ibface]
end


"""
    $(SIGNATURES)

Fill edge with the help of assemblydata.
"""
function _fill!(
    edge::Edge,
    asmdata::CellwiseAssemblyData{Tv,Ti},
    iedge,
    icell,
) where {Tv,Ti}
    if edge.has_celledges #  cellx==celledges, edgenodes==global_edgenodes
        # If we work with projections of fluxes onto edges,
        # we need to ensure that the edges are accessed with the
        # same orientation without regard of the orientation induced
        # by local cell numbering
        edge.index = edge.cellx[iedge, icell]
        edge.node[1] = edge.edgenodes[1, edge.index]
        edge.node[2] = edge.edgenodes[2, edge.index]
    else # cx==cellnodes, edgenodes== local_edgenodes
        edge.index = 0
        edge.node[1] = edge.cellx[edge.edgenodes[1, iedge], icell]
        edge.node[2] = edge.cellx[edge.edgenodes[2, iedge], icell]
    end
    edge.region = edge.cellregions[icell]
    edge.fac = asmdata.edgefactors[iedge, icell]
    edge.icell = icell
end



"""
    $(SIGNATURES)

Fill boundary edge with the help of assemblydata.
"""
function _fill!(
    bedge::BEdge,
    asmdata::CellwiseAssemblyData{Tv,Ti},
    ibedge,
    ibface,
) where {Tv,Ti}
    bedge.index = bedge.bfaceedges[ibedge, ibface]
    bedge.node[1] = bedge.bedgenodes[1, bedge.index]
    bedge.node[2] = bedge.bedgenodes[2, bedge.index]
    bedge.region = bedge.bfaceregions[ibface]
    bedge.icell = ibface
    bedge.fac = asmdata.edgefactors[ibedge, ibface]
end


"""
    $(SIGNATURES)

Fill node with the help of assemblydata.
"""
function _fill!(node::Node, asmdata::EdgewiseAssemblyData{Tv,Ti}, k, inode) where {Tv,Ti}
    node.index = inode
    node.region = asmdata.nodefactors.rowval[k]
    node.fac = asmdata.nodefactors.nzval[k]
end

"""
    $(SIGNATURES)

Fill edge with the help of assemblydata.
"""
function _fill!(edge::Edge, asmdata::EdgewiseAssemblyData{Tv,Ti}, k, iedge) where {Tv,Ti}
    edge.index = iedge
    edge.node[1] = edge.edgenodes[1, edge.index]
    edge.node[2] = edge.edgenodes[2, edge.index]
    edge.region = asmdata.edgefactors.rowval[k]
    edge.fac = asmdata.edgefactors.nzval[k]
end





"""
$(SIGNATURES)

Assemble residual and jacobian for node functions. Parameters:

- `system`: System to be worked with
- `node`: node

- `asm_jac(idof,jdof,ispec,jspec)`: e.g.  assemble entry `ispec,jspec` of local jacobian into entry `idof,jdof` of global matrix
- `asm_param(idof,ispec,iparam)` shall assemble parameter derivatives
"""
@inline function assemble_res_jac(
    node::Node,
    system::AbstractSystem,
    asm_res::R,
    asm_jac::J,
    asm_param::P,
) where {R,J,P}
    K = node.index
    ireg = node.region
    for idof = firstnodedof(system, K):lastnodedof(system, K)
        ispec = getspecies(system, idof)
        if isregionspecies(system, ispec, ireg) # it is not enough to know if the species are defined...
            asm_res(idof, ispec)
            for jdof = firstnodedof(system, K):lastnodedof(system, K)
                jspec = getspecies(system, jdof)
                if isregionspecies(system, jspec, ireg)
                    asm_jac(idof, jdof, ispec, jspec)
                end
            end
            for iparam = 1:(system.num_parameters)
                asm_param(idof, ispec, iparam)
            end
        end
    end
end

"""
$(SIGNATURES)

Assemble residual and jacobian for boundary node functions.
See [`assemble_res_jac`](@ref) for more explanations.
"""
@inline function assemble_res_jac(
    bnode::BNode,
    system::AbstractSystem,
    asm_res::R,
    asm_jac::J,
    asm_param::P,
) where {R,J,P}
    K = bnode.index
    for idof = firstnodedof(system, K):lastnodedof(system, K)
        ispec = getspecies(system, idof)
        if isnodespecies(system, ispec, K)
            asm_res(idof, ispec)
            for jdof = firstnodedof(system, K):lastnodedof(system, K)
                jspec = getspecies(system, jdof)
                if isnodespecies(system, jspec, K)
                    asm_jac(idof, jdof, ispec, jspec)
                end
                for iparam = 1:(system.num_parameters)
                    asm_param(idof, ispec, iparam)
                end
            end
        end
    end
end

"""
$(SIGNATURES)

Assemble residual for node functions.
See [`assemble_res_jac`](@ref) for more explanations.
"""
@inline function assemble_res(node::Node, system::AbstractSystem, asm_res::R) where {R}
    K = node.index
    ireg = node.region
    for idof = firstnodedof(system, K):lastnodedof(system, K)
        ispec = getspecies(system, idof)
        if isregionspecies(system, ispec, ireg)
            asm_res(idof, ispec)
        end
    end
end

"""
$(SIGNATURES)

Assemble residual for boundary node functions.
See [`assemble_res_jac`](@ref) for more explanations.
"""
@inline function assemble_res(bnode::BNode, system::AbstractSystem, asm_res::R) where {R}
    K = bnode.index
    for idof = firstnodedof(system, K):lastnodedof(system, K)
        ispec = getspecies(system, idof)
        if isnodespecies(system, ispec, K)
            asm_res(idof, ispec)
        end
    end
end

"""
$(SIGNATURES)

Assemble residual and jacobian for edge (flux) functions. Parameters:

- `system`: System to be worked with
- `edge`: edge
- `asm_res(idofK,idofL,ispec)`: e.g. assemble local ispec to global degrees of freedom in unknowns
- `asm_jac(idofK,jdofK,idofL,jdofL,ispec,jspec)`: e.g.  assemble entry `ispec,jspec` of local jacobian into entry four entries defined by `idofK` and `idofL` of global matrix
- `asm_param(idofK,idofL,ispec,iparam)` shall assemble parameter derivatives
"""
@inline function assemble_res_jac(
    edge::Edge,
    system::AbstractSystem,
    asm_res::R,
    asm_jac::J,
    asm_param::P,
) where {R,J,P}
    K = edge.node[1]
    L = edge.node[2]
    ireg = edge.region
    for idofK = firstnodedof(system, K):lastnodedof(system, K)
        ispec = getspecies(system, idofK)
        if isregionspecies(system, ispec, ireg)
            idofL = getnodedof(system, ispec, L)
            if idofL > 0
                asm_res(idofK, idofL, ispec)
                for jdofK = firstnodedof(system, K):lastnodedof(system, K)
                    jspec = getspecies(system, jdofK)
                    if isregionspecies(system, jspec, ireg)
                        jdofL = getnodedof(system, jspec, L)
                        if jdofL > 0
                            asm_jac(idofK, jdofK, idofL, jdofL, ispec, jspec)
                        end
                    end
                end
                for iparam = 1:(system.num_parameters)
                    asm_param(idofK, idofL, ispec, iparam)
                end
            end
        end

    end
end

"""
$(SIGNATURES)

Assemble residual for edge (flux) functions.
See [`assemble_res_jac`](@ref) for more explanations.
"""
@inline function assemble_res(edge::Edge, system::AbstractSystem, asm_res::R) where {R}
    K = edge.node[1]
    L = edge.node[2]
    ireg = edge.region

    for idofK = firstnodedof(system, K):lastnodedof(system, K)
        ispec = getspecies(system, idofK)
        if isregionspecies(system, ispec, ireg)
            idofL = getnodedof(system, ispec, L)
            if idofL > 0
                asm_res(idofK, idofL, ispec)
            end
        end
    end
end

"""
$(SIGNATURES)

Assemble residual and jacobian for boundary edge (flux) functions.
See [`assemble_res_jac`](@ref) for more explanations.
"""
@inline function assemble_res_jac(
    bedge::BEdge,
    system::AbstractSystem,
    asm_res::R,
    asm_jac::J,
    asm_param::P,
) where {R,J,P}
    K = bedge.node[1]
    L = bedge.node[2]
    for idofK = firstnodedof(system, K):lastnodedof(system, K)
        ispec = getspecies(system, idofK)
        if isnodespecies(system, ispec, K)
            idofL = getnodedof(system, ispec, L)
            if idofL > 0
                asm_res(idofK, idofL, ispec)

                for jdofK = firstnodedof(system, K):lastnodedof(system, K)
                    jspec = getspecies(system, jdofK)
                    if isnodespecies(system, jspec, K)
                        jdofL = getnodedof(system, jspec, L)
                        if jdofL > 0
                            asm_jac(idofK, jdofK, idofL, jdofL, ispec, jspec)
                        end
                    end
                end
                for iparam = 1:(system.num_parameters)
                    asm_param(idofK, idofL, ispec, iparam)
                end
            end
        end
    end

end

"""
$(SIGNATURES)

Assemble residual for boundary edge (flux) functions.
See [`assemble_res_jac`](@ref) for more explanations.
"""
@inline function assemble_res(bedge::BEdge, system::AbstractSystem, asm_res::R) where {R}
    K = bedge.node[1]
    L = bedge.node[2]
    for idofK = firstnodedof(system, K):lastnodedof(system, K)
        ispec = getspecies(system, idofK)
        if isnodespecies(system, ispec, K)
            idofL = dof(F, ispec, L)
            if idofL > 0
                asm_res(idofK, idofL, ispec)
            end
        end
    end
end
