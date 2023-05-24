abstract type AbstractAssemblyData{Tv,Ti} end


struct CellWiseAssemblyData{Tv,Ti} <: AbstractAssemblyData{Tv,Ti}
    """
    Precomputed geometry factors for cell nodes
    """
    nodefactors::Array{Tv,2}

    """
    Precomputed geometry factors for cell edges
    """
    edgefactors::Array{Tv,2}

    """
    Cell regions
    """
    cellregions::Array{Ti,1}
end


struct EdgeWiseAssemblyData{Tv,Ti} <: AbstractAssemblyData{Tv,Ti}
    """
        Precomputed geometry factors for cell nodes
    """
    nodefactors::SparseMatrixCSC{Tv,Ti}

    """
    Precomputed geometry factors for cell edges
    """
    edgefactors::SparseMatrixCSC{Tv,Ti}

end


nodebatch(asmdata::CellWiseAssemblyData{Tv,Ti}) where {Tv,Ti} =
    1:size(asmdata.nodefactors, 2)
edgebatch(asmdata::CellWiseAssemblyData{Tv,Ti}) where {Tv,Ti} =
    1:size(asmdata.edgefactors, 2)

nodebatch(asmdata::EdgeWiseAssemblyData{Tv,Ti}) where {Tv,Ti} =
    1:size(asmdata.nodefactors, 2)
edgebatch(asmdata::EdgeWiseAssemblyData{Tv,Ti}) where {Tv,Ti} =
    1:size(asmdata.edgefactors, 2)


noderange(asmdata::EdgeWiseAssemblyData{Tv,Ti}, inode) where {Tv,Ti} =
    nzrange(asmdata.nodefactors, inode)
edgerange(asmdata::EdgeWiseAssemblyData{Tv,Ti}, iedge) where {Tv,Ti} =
    nzrange(asmdata.edgefactors, iedge)

noderange(asmdata::CellWiseAssemblyData{Tv,Ti}, inode) where {Tv,Ti} =
    1:size(asmdata.nodefactors, 1)
edgerange(asmdata::CellWiseAssemblyData{Tv,Ti}, iedge) where {Tv,Ti} =
    1:size(asmdata.edgefactors, 1)


function _fill!(
    node::Node,
    asmdata::CellWiseAssemblyData{Tv,Ti},
    inode,
    icell,
) where {Tv,Ti}
    node.fac = asmdata.nodefactors[inode, icell]
    _fill!(node, inode, icell)
end

function _fill!(
    edge::Edge,
    asmdata::CellWiseAssemblyData{Tv,Ti},
    iedge,
    icell,
) where {Tv,Ti}
    edge.fac = asmdata.edgefactors[iedge, icell]
    _fill!(edge, iedge, icell)
end





function _fill!(node::Node, asmdata::EdgeWiseAssemblyData{Tv,Ti}, k, inode) where {Tv,Ti}
    _xfill!(node, asmdata.nodefactors.rowval[k], inode)
    node.fac = asmdata.nodefactors.nzval[k]
end

function _fill!(edge::Edge, asmdata::EdgeWiseAssemblyData{Tv,Ti}, k, iedge) where {Tv,Ti}
    _xfill!(edge, asmdata.edgefactors.rowval[k], iedge)
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
        end
        for iparam = 1:(system.num_parameters)
            asm_param(idof, ispec, iparam)
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
            end
        end
        for iparam = 1:(system.num_parameters)
            asm_param(idof, ispec, iparam)
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
- `node`: node
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
            end
        end

        for iparam = 1:(system.num_parameters)
            asm_param(idofK, idofL, ispec, iparam)
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
            end
        end
    end

    for iparam = 1:(system.num_parameters)
        asm_param(idofK, idofL, ispec, iparam)
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
