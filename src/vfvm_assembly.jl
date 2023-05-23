##################################################################
"""
$(SIGNATURES)

Extract value from dual number. Use to debug physics callbacks.
Re-exported from ForwardDiff.jl
"""
const value = ForwardDiff.value


"""
$(SIGNATURES)

Add value `v*fac` to matrix if `v` is nonzero
"""
@inline function _addnz(matrix, i, j, v::Tv, fac) where {Tv}
    if isnan(v)
        error("trying to assemble NaN")
    end
    if v != zero(Tv)
        rawupdateindex!(matrix, +, v * fac, i, j)
    end
end

ExtendableSparse.rawupdateindex!(m::AbstractMatrix, op, v, i, j) = m[i, j] = op(m[i, j], v)


function zero!(m::ExtendableSparseMatrix{Tv,Ti}) where {Tv,Ti}
    nzv = nonzeros(m)
    nzv .= zero(Tv)
end

zero!(m::AbstractMatrix{T}) where {T} = m .= zero(T)

"""
$(SIGNATURES)

Main assembly method.

Evaluate solution with result in right hand side F and 
assemble Jacobi matrix into system.matrix.
"""
function eval_and_assemble(
    system::System{Tv,Tc,Ti,Tm,TSpecMat,TSolArray},
    U::AbstractMatrix{Tv}, # Actual solution iteration
    UOld::AbstractMatrix{Tv}, # Old timestep solution
    F::AbstractMatrix{Tv},# Right hand side
    time,
    tstep,# time step size. Inf means stationary solution
    λ,
    params::AbstractVector;
    edge_cutoff = 0.0,
) where {Tv,Tc,Ti,Tm,TSpecMat,TSolArray}
    _complete!(system) # needed here as well for test function system which does not use newton

    grid = system.grid
    physics = system.physics
    node = Node(system, time, λ, params)
    edge = Edge(system, time, λ, params)
    nspecies::Int = num_species(system)

    # Reset matrix + rhs
    zero!(system.matrix)
    F .= 0.0
    nparams::Int = system.num_parameters

    dudp = system.dudp

    for iparam = 1:nparams
        dudp[iparam] .= 0.0
    end

    # Arrays for gathering solution data
    UK = Array{Tv,1}(undef, nspecies + nparams)
    UKOld = Array{Tv,1}(undef, nspecies + nparams)
    UKL = Array{Tv,1}(undef, 2 * nspecies + nparams)

    @assert length(params) == nparams
    if nparams > 0
        UK[(nspecies+1):end] .= params
        UKOld[(nspecies+1):end] .= params
        UKL[(2*nspecies+1):end] .= params
    end

    # Inverse of timestep
    # According to Julia documentation, 1/Inf=0 which
    # comes handy to write compact code here for the
    # case of stationary problems.
    tstepinv = 1.0 / tstep

    #
    # These wrap the different physics functions.
    #
    src_evaluator = ResEvaluator(physics, :source, UK, node, nspecies)
    rea_evaluator = ResJacEvaluator(physics, :reaction, UK, node, nspecies)
    stor_evaluator = ResJacEvaluator(physics, :storage, UK, node, nspecies)
    oldstor_evaluator = ResEvaluator(physics, :storage, UK, node, nspecies)
    flux_evaluator = ResJacEvaluator(physics, :flux, UKL, edge, nspecies)
    erea_evaluator = ResJacEvaluator(physics, :edgereaction, UKL, edge, nspecies)


    ncalloc = @allocated for icell in noderange(system.assembly_data)
        for inode in noderange(system.assembly_data,icell)
            _fill!(node,system.assembly_data,inode,icell)
            fac=node.fac
                        @views UK[1:nspecies] .= U[:, node.index]
            @views UKOld[1:nspecies] .= UOld[:, node.index]

            evaluate!(src_evaluator)
            src = res(src_evaluator)

            evaluate!(rea_evaluator, UK)
            res_react = res(rea_evaluator)
            jac_react = jac(rea_evaluator)

            evaluate!(stor_evaluator, UK)
            res_stor = res(stor_evaluator)
            jac_stor = jac(stor_evaluator)

            evaluate!(oldstor_evaluator, UKOld)
            oldstor = res(oldstor_evaluator)

            @inline function asm_res(idof, ispec)
                _add(
                    F,
                    idof,
                    fac * (
                        res_react[ispec] - src[ispec] +
                        (res_stor[ispec] - oldstor[ispec]) * tstepinv
                    ),
                )
            end

            @inline function asm_jac(idof, jdof, ispec, jspec)
                _addnz(
                    system.matrix,
                    idof,
                    jdof,
                    jac_react[ispec, jspec] + jac_stor[ispec, jspec] * tstepinv,
                    fac,
                )
            end

            @inline function asm_param(idof, ispec, iparam)
                jparam = nspecies + iparam
                dudp[iparam][ispec, idof] +=
                    (jac_react[ispec, jparam] + jac_stor[ispec, jparam] * tstepinv) * fac
            end

            assemble_res_jac(node, system, asm_res, asm_jac, asm_param)

        end
    end
    
    ncalloc += @allocated for icell in edgerange(system.assembly_data)
        for iedge in edgerange(system.assembly_data,icell)
            _fill!(edge,system.assembly_data,iedge,icell) 
            fac=edge.fac
            @views UKL[1:nspecies] .= U[:, edge.node[1]]
            @views UKL[(nspecies+1):(2*nspecies)] .= U[:, edge.node[2]]
            
            evaluate!(flux_evaluator, UKL)
            res_flux = res(flux_evaluator)
            jac_flux = jac(flux_evaluator)
            
            @inline function asm_res(idofK, idofL, ispec)
                val = fac * res_flux[ispec]
                _add(F, idofK, val)
                _add(F, idofL, -val)
            end
            
            @inline function asm_jac(idofK, jdofK, idofL, jdofL, ispec, jspec)
                _addnz(system.matrix, idofK, jdofK, +jac_flux[ispec, jspec], fac)
                _addnz(system.matrix, idofL, jdofK, -jac_flux[ispec, jspec], fac)
                _addnz(system.matrix, idofK, jdofL, +jac_flux[ispec, jspec+nspecies], fac)
                _addnz(system.matrix, idofL, jdofL, -jac_flux[ispec, jspec+nspecies], fac)
            end
            
            @inline function asm_param(idofK, idofL, ispec, iparam)
                jparam = 2 * nspecies + iparam
                dudp[iparam][ispec, idofK] += fac * jac_flux[ispec, jparam]
                dudp[iparam][ispec, idofL] -= fac * jac_flux[ispec, jparam]
            end
            
            assemble_res_jac(edge, system, asm_res, asm_jac, asm_param)
            
            if isnontrivial(erea_evaluator)
                evaluate!(erea_evaluator, UKL)
                res_erea = res(erea_evaluator)
                jac_erea = jac(erea_evaluator)
                
                @inline function ereaasm_res(idofK, idofL, ispec)
                    val = fac * res_erea[ispec]
                    _add(F, idofK, val)
                    _add(F, idofL, val)
                end
                
                @inline function ereaasm_jac(idofK, jdofK, idofL, jdofL, ispec, jspec)
                    _addnz(system.matrix, idofK, jdofK, +jac_erea[ispec, jspec], fac)
                    _addnz(system.matrix, idofL, jdofK, -jac_erea[ispec, jspec], fac)
                    _addnz(
                        system.matrix,
                        idofK,
                        jdofL,
                        -jac_erea[ispec, jspec+nspecies],
                        fac,
                    )
                    _addnz(
                        system.matrix,
                        idofL,
                        jdofL,
                        +jac_erea[ispec, jspec+nspecies],
                        fac,
                    )
                end
                
                @inline function ereaasm_param(idofK, idofL, ispec, iparam)
                    jparam = 2 * nspecies + iparam
                    dudp[iparam][ispec, idofK] += fac * jac_erea[ispec, jparam]
                    dudp[iparam][ispec, idofL] += fac * jac_erea[ispec, jparam]
                end
                
                assemble_res_jac(edge, system, ereaasm_res, ereaasm_jac, ereaasm_param)
            end
        end
    end

    bnode = BNode(system, time, λ, params)
    bedge = BEdge(system, time, λ, params)


    bfacenodefactors::Array{Tv,2} = system.bfacenodefactors
    bfaceedgefactors::Array{Tv,2} = system.bfaceedgefactors
    boundary_factors::Array{Tv,2} = system.boundary_factors
    boundary_values::Array{Tv,2} = system.boundary_values
    has_legacy_bc = !iszero(boundary_factors) || !iszero(boundary_values)

    bsrc_evaluator = ResEvaluator(physics, :bsource, UK, bnode, nspecies)
    brea_evaluator = ResJacEvaluator(physics, :breaction, UK, bnode, nspecies)
    bstor_evaluator = ResJacEvaluator(physics, :bstorage, UK, bnode, nspecies)
    oldbstor_evaluator = ResEvaluator(physics, :bstorage, UK, bnode, nspecies)
    bflux_evaluator = ResJacEvaluator(physics, :bflux, UKL, bedge, nspecies)

    
    bfaceregions::Vector{Ti} = grid[BFaceRegions]
    nbfaces = num_bfaces(grid)
    bgeom = grid[BFaceGeometries][1]

    nbn::Int = num_nodes(bgeom)
    nbe::Int = num_edges(bgeom)

    # Assembly loop for boundary conditions
    nballoc = @allocated for ibface = 1:nbfaces

        ibreg = bfaceregions[ibface]

        # Loop over nodes of boundary face
        for ibnode = 1:nbn
            # Fill bnode data shuttle with data from grid
            _fill!(bnode, ibnode, ibface)

            # Measure of boundary face part assembled to node
            bnode_factor::Tv = bfacenodefactors[ibnode, ibface]

            if has_legacy_bc
                # Global index of node
                K = bnode.index
                bnode.Dirichlet = Dirichlet / bnode_factor
                # Assemble "standard" boundary conditions: Robin or
                # Dirichlet
                # valid only for interior species, currently not checked
                for ispec = 1:nspecies
                    idof = dof(F, ispec, K)
                    # If species is present, assemble the boundary condition
                    if idof > 0
                        # Get user specified data
                        boundary_factor = boundary_factors[ispec, ibreg]
                        boundary_value = boundary_values[ispec, ibreg]

                        if boundary_factor == Dirichlet
                            # Dirichlet is encoded in the boundary condition factor
                            # Penalty method: scale   (u-boundary_value) by penalty=boundary_factor

                            # Add penalty*boundary_value to right hand side
                            F[ispec, K] += boundary_factor * (U[ispec, K] - boundary_value)

                            # Add penalty to matrix main diagonal (without bnode factor, so penalty
                            # is independent of h)
                            _addnz(system.matrix, idof, idof, boundary_factor, 1)
                        else
                            # Robin boundary condition
                            F[ispec, K] +=
                                bnode_factor *
                                (boundary_factor * U[ispec, K] - boundary_value)
                            _addnz(system.matrix, idof, idof, boundary_factor, bnode_factor)
                        end
                    end
                end
            end

            # Copy unknown values from solution into dense array
            @views UK[1:nspecies] .= U[:, bnode.index]

            evaluate!(bsrc_evaluator)
            bsrc = res(bsrc_evaluator)

            evaluate!(brea_evaluator, UK)
            res_breact = res(brea_evaluator)
            jac_breact = jac(brea_evaluator)

            asm_res1(idof, ispec) =
                _add(F, idof, bnode_factor * (res_breact[ispec] - bsrc[ispec]))

            asm_jac1(idof, jdof, ispec, jspec) = _addnz(
                system.matrix,
                idof,
                jdof,
                jac_breact[ispec, jspec],
                bnode_factor,
            )

            asm_param1(idof, ispec, iparam) =
                dudp[iparam][ispec, idof] +=
                    jac_breact[ispec, nspecies+iparam] * bnode_factor

            assemble_res_jac(bnode, system, asm_res1, asm_jac1, asm_param1)

            if isnontrivial(bstor_evaluator)
                evaluate!(bstor_evaluator, UK)
                res_bstor = res(bstor_evaluator)
                jac_bstor = jac(bstor_evaluator)

                @views UKOld .= UOld[:, bnode.index]
                evaluate!(oldbstor_evaluator, UKOld)
                oldbstor = res(oldbstor_evaluator)

                asm_res2(idof, ispec) = _add(
                    F,
                    idof,
                    bnode_factor * (res_bstor[ispec] - oldbstor[ispec]) * tstepinv,
                )

                function asm_jac2(idof, jdof, ispec, jspec)
                    _addnz(
                        system.matrix,
                        idof,
                        jdof,
                        jac_bstor[ispec, jspec],
                        bnode_factor * tstepinv,
                    )
                end

                function asm_param2(idof, ispec, iparam)
                    dudp[iparam][ispec, idof] +=
                        jac_bstor[ispec, nspecies+iparam] * bnode_factor * tstepinv
                end

                assemble_res_jac(bnode, system, asm_res2, asm_jac2, asm_param2)
            end
        end # ibnode=1:nbn

        if isnontrivial(bflux_evaluator)
            for ibedge = 1:nbe
                fac = bfaceedgefactors[ibedge, ibface]
                #if abs(fac) < edge_cutoff    continue    end

                _fill!(bedge, ibedge, ibface)
                @views UKL[1:nspecies] .= U[:, bedge.node[1]]
                @views UKL[(nspecies+1):(2*nspecies)] .= U[:, bedge.node[2]]

                evaluate!(bflux_evaluator, UKL)
                res_bflux = res(bflux_evaluator)
                jac_bflux = jac(bflux_evaluator)

                function asm_res(idofK, idofL, ispec)
                    _add(F, idofK, fac * res_bflux[ispec])
                    _add(F, idofL, -fac * res_bflux[ispec])
                end

                function asm_jac(idofK, jdofK, idofL, jdofL, ispec, jspec)
                    _addnz(system.matrix, idofK, jdofK, +jac_bflux[ispec, jspec], fac)
                    _addnz(system.matrix, idofL, jdofK, -jac_bflux[ispec, jspec], fac)
                    _addnz(
                        system.matrix,
                        idofK,
                        jdofL,
                        +jac_bflux[ispec, jspec+nspecies],
                        fac,
                    )
                    _addnz(
                        system.matrix,
                        idofL,
                        jdofL,
                        -jac_bflux[ispec, jspec+nspecies],
                        fac,
                    )
                end

                function asm_param(idofK, idofL, ispec, iparam)
                    jparam = 2 * nspecies + iparam
                    dudp[iparam][ispec, idofK] += fac * jac_bflux[ispec, jparam]
                    dudp[iparam][ispec, idofL] -= fac * jac_bflux[ispec, jparam]
                end
                assemble_res_jac(bedge, system, asm_res, asm_jac, asm_param)
            end
        end
    end
    noallocs(m::ExtendableSparseMatrix) = isnothing(m.lnkmatrix)
    noallocs(m::AbstractMatrix) = false
    # if  no new matrix entries have been created, we should see no allocations
    # in the previous two loops
    neval = 1
    if !noallocs(system.matrix)
        ncalloc = 0
        nballoc = 0
        neval = 0
    end
    _eval_and_assemble_generic_operator(system, U, F)
    _eval_and_assemble_inactive_species(system, U, UOld, F)

    ncalloc, nballoc, neval
end

"""
Evaluate and assemble jacobian for generic operator part.
"""
function _eval_and_assemble_generic_operator(system::AbstractSystem, U, F)
    if !has_generic_operator(system)
        return
    end
    generic_operator(f, u) = system.physics.generic_operator(f, u, system)
    vecF = values(F)
    vecU = values(U)
    y = similar(vecF)
    generic_operator(y, vecU)
    vecF .+= y
    forwarddiff_color_jacobian!(
        system.generic_matrix,
        generic_operator,
        vecU;
        colorvec = system.generic_matrix_colors,
    )
    rowval = system.generic_matrix.rowval
    colptr = system.generic_matrix.colptr
    nzval = system.generic_matrix.nzval
    for i = 1:(length(colptr)-1)
        for j = colptr[i]:(colptr[i+1]-1)
            updateindex!(system.matrix, +, nzval[j], rowval[j], i)
        end
    end
end

function _eval_generic_operator(system::AbstractSystem, U, F)
    if !has_generic_operator(system)
        return
    end
    generic_operator(f, u) = system.physics.generic_operator(f, u, system)
    vecF = values(F)
    vecU = values(U)
    y = similar(vecF)
    generic_operator(y, vecU)
    vecF .+= y
end
