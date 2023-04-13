##################################################################
"""
$(SIGNATURES)

Extract value from dual number. Use to debug physics callbacks.
Re-exported from ForwardDiff.jl
"""
const value = ForwardDiff.value

# These are needed to enable iterative solvers to work with dual numbers
Base.Float64(x::ForwardDiff.Dual) = value(x)
function Random.rand(rng::AbstractRNG,
                     ::Random.SamplerType{ForwardDiff.Dual{T, V, N}}) where {T, V, N}
    ForwardDiff.Dual{T, V, N}(rand(rng, V))
end

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

"""
$(TYPEDEF)

Exception thrown if Newton's method convergence fails.
"""
struct ConvergenceError <: Exception end

"""
$(TYPEDEF)

Exception thrown if error occured during assembly (e.g. domain error)
"""
struct AssemblyError <: Exception end

"""
$(TYPEDEF)

Exception thrown if error occured during factorization.
"""
struct LinearSolverError <: Exception end

"""
$(TYPEDEF)

Exception thrown if embedding fails
"""
struct EmbeddingError <: Exception
    msg::String
end

"""
Print error when catching exceptions
"""
function _print_error(err, st)
    println()
    println(err)
    nlines = 5
    for i = 1:min(nlines, length(st))
        line = @sprintf("%s", st[i])
        L = length(line)
        if L < 80
            println(line)
        else
            print(line[1:35])
            print(" ... ")
            println(line[(L - 35):L])
        end
    end
    if length(st) > nlines
        println("...")
    end
    println()
end

function zero!(m::ExtendableSparseMatrix{Tv, Ti}) where {Tv, Ti}
    nzv = nonzeros(m)
    nzv .= zero(Tv)
end

zero!(m::AbstractMatrix{T}) where {T} = m .= zero(T)

################################################################
"""
$(SIGNATURES)

Main assembly method.

Evaluate solution with result in right hand side F and 
assemble Jacobi matrix into system.matrix.
"""
function eval_and_assemble(system::System{Tv, Tc, Ti, Tm, TSpecMat, TSolArray},
                           U::AbstractMatrix{Tv}, # Actual solution iteration
                           UOld::AbstractMatrix{Tv}, # Old timestep solution
                           F::AbstractMatrix{Tv},# Right hand side
                           time,
                           tstep,# time step size. Inf means stationary solution
                           λ,
                           params::AbstractVector;
                           edge_cutoff = 0.0) where {Tv, Tc, Ti, Tm, TSpecMat, TSolArray}
    _complete!(system) # needed here as well for test function system which does not use newton

    grid = system.grid
    physics = system.physics
    node = Node(system, time, λ, params)
    bnode = BNode(system, time, λ, params)
    edge = Edge(system, time, λ, params)
    bedge = BEdge(system, time, λ, params)

    nspecies::Int = num_species(system)

    cellnodefactors::Array{Tv, 2} = system.cellnodefactors
    celledgefactors::Array{Tv, 2} = system.celledgefactors
    bfacenodefactors::Array{Tv, 2} = system.bfacenodefactors
    bfaceedgefactors::Array{Tv, 2} = system.bfaceedgefactors

    # Reset matrix + rhs
    zero!(system.matrix)
    F .= 0.0
    nparams::Int = system.num_parameters

    dudp = system.dudp

    for iparam = 1:nparams
        dudp[iparam] .= 0.0
    end

    # Arrays for gathering solution data
    UK = Array{Tv, 1}(undef, nspecies + nparams)
    UKOld = Array{Tv, 1}(undef, nspecies + nparams)
    UKL = Array{Tv, 1}(undef, 2 * nspecies + nparams)

    @assert length(params) == nparams
    if nparams > 0
        UK[(nspecies + 1):end] .= params
        UKOld[(nspecies + 1):end] .= params
        UKL[(2 * nspecies + 1):end] .= params
    end

    # Inverse of timestep
    # According to Julia documentation, 1/Inf=0 which
    # comes handy to write compact code here for the
    # case of stationary problems.
    tstepinv = 1.0 / tstep

    boundary_factors::Array{Tv, 2} = system.boundary_factors
    boundary_values::Array{Tv, 2} = system.boundary_values
    bfaceregions::Vector{Ti} = grid[BFaceRegions]
    cellregions::Vector{Ti} = grid[CellRegions]

    nbfaces = num_bfaces(grid)
    ncells = num_cells(grid)
    geom = grid[CellGeometries][1]
    bgeom = grid[BFaceGeometries][1]

    nn::Int = num_nodes(geom)
    ne::Int = num_edges(geom)

    has_legacy_bc = !iszero(boundary_factors) || !iszero(boundary_values)

    #
    # These wrap the different physics functions.
    #
    src_evaluator = ResEvaluator(physics, :source, UK, node, nspecies)
    rea_evaluator = ResJacEvaluator(physics, :reaction, UK, node, nspecies)
    stor_evaluator = ResJacEvaluator(physics, :storage, UK, node, nspecies)
    oldstor_evaluator = ResEvaluator(physics, :storage, UK, node, nspecies)
    flux_evaluator = ResJacEvaluator(physics, :flux, UKL, edge, nspecies)
    erea_evaluator = ResJacEvaluator(physics, :edgereaction, UKL, edge, nspecies)

    bsrc_evaluator = ResEvaluator(physics, :bsource, UK, bnode, nspecies)
    brea_evaluator = ResJacEvaluator(physics, :breaction, UK, bnode, nspecies)
    bstor_evaluator = ResJacEvaluator(physics, :bstorage, UK, bnode, nspecies)
    oldbstor_evaluator = ResEvaluator(physics, :bstorage, UK, bnode, nspecies)
    bflux_evaluator = ResJacEvaluator(physics, :bflux, UKL, bedge, nspecies)

    ncalloc = @allocated for icell = 1:ncells
        ireg = cellregions[icell]
        for inode = 1:nn
            _fill!(node, inode, icell)
            @views UK[1:nspecies] .= U[:, node.index]
            @views UKOld[1:nspecies] .= UOld[:, node.index]
            fac = cellnodefactors[inode, icell]

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
                _add(F,
                     idof,
                     fac * (res_react[ispec] - src[ispec] +
                            (res_stor[ispec] - oldstor[ispec]) * tstepinv))
            end

            @inline function asm_jac(idof, jdof, ispec, jspec)
                _addnz(system.matrix,
                       idof,
                       jdof,
                       jac_react[ispec, jspec] + jac_stor[ispec, jspec] * tstepinv,
                       fac)
            end

            @inline function asm_param(idof, ispec, iparam)
                jparam = nspecies + iparam
                dudp[iparam][ispec, idof] += (jac_react[ispec, jparam] + jac_stor[ispec, jparam] * tstepinv) * fac
            end

            assemble_res_jac(node, system, asm_res, asm_jac, asm_param)
        end

        for iedge = 1:ne
            abs(celledgefactors[iedge, icell]) < edge_cutoff && continue
            _fill!(edge, iedge, icell)
            fac = celledgefactors[iedge, icell]

            #Set up argument for fluxwrap
            @views UKL[1:nspecies] .= U[:, edge.node[1]]
            @views UKL[(nspecies + 1):(2 * nspecies)] .= U[:, edge.node[2]]

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
                _addnz(system.matrix, idofK, jdofL, +jac_flux[ispec, jspec + nspecies], fac)
                _addnz(system.matrix, idofL, jdofL, -jac_flux[ispec, jspec + nspecies], fac)
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
                    _addnz(system.matrix, idofK, jdofL, -jac_erea[ispec, jspec + nspecies], fac)
                    _addnz(system.matrix, idofL, jdofL, +jac_erea[ispec, jspec + nspecies], fac)
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

                            # Add penalty*boundry_value to right hand side
                            F[ispec, K] += boundary_factor * (U[ispec, K] - boundary_value)

                            # Add penalty to matrix main diagonal (without bnode factor, so penalty
                            # is independent of h)
                            _addnz(system.matrix, idof, idof, boundary_factor, 1)
                        else
                            # Robin boundary condition
                            F[ispec, K] += bnode_factor *
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

            asm_res1(idof, ispec) = _add(F, idof, bnode_factor * (res_breact[ispec] - bsrc[ispec]))

            asm_jac1(idof, jdof, ispec, jspec) = _addnz(system.matrix,
                                                        idof,
                                                        jdof,
                                                        jac_breact[ispec, jspec],
                                                        bnode_factor)

            asm_param1(idof, ispec, iparam) = dudp[iparam][ispec, idof] += jac_breact[ispec, nspecies + iparam] * bnode_factor

            assemble_res_jac(bnode, system, asm_res1, asm_jac1, asm_param1)

            if isnontrivial(bstor_evaluator)
                evaluate!(bstor_evaluator, UK)
                res_bstor = res(bstor_evaluator)
                jac_bstor = jac(bstor_evaluator)

                @views UKOld .= UOld[:, bnode.index]
                evaluate!(oldbstor_evaluator, UKOld)
                oldbstor = res(oldbstor_evaluator)

                asm_res2(idof, ispec) = _add(F,
                                             idof,
                                             bnode_factor * (res_bstor[ispec] - oldbstor[ispec]) * tstepinv)

                function asm_jac2(idof, jdof, ispec, jspec)
                    _addnz(system.matrix,
                           idof,
                           jdof,
                           jac_bstor[ispec, jspec],
                           bnode_factor * tstepinv)
                end

                function asm_param2(idof, ispec, iparam)
                    dudp[iparam][ispec, idof] += jac_bstor[ispec, nspecies + iparam] * bnode_factor * tstepinv
                end

                assemble_res_jac(bnode, system, asm_res2, asm_jac2, asm_param2)
            end
        end # ibnode=1:nbn

        if isnontrivial(bflux_evaluator)
            for ibedge = 1:nbe
                fac = bfaceedgefactors[ibedge, ibface]
                if abs(fac) < edge_cutoff
                    continue
                end

                _fill!(bedge, ibedge, ibface)
                @views UKL[1:nspecies] .= U[:, bedge.node[1]]
                @views UKL[(nspecies + 1):(2 * nspecies)] .= U[:, bedge.node[2]]

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
                    _addnz(system.matrix,
                           idofK,
                           jdofL,
                           +jac_bflux[ispec, jspec + nspecies],
                           fac)
                    _addnz(system.matrix,
                           idofL,
                           jdofL,
                           -jac_bflux[ispec, jspec + nspecies],
                           fac)
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
    forwarddiff_color_jacobian!(system.generic_matrix,
                                generic_operator,
                                vecU;
                                colorvec = system.generic_matrix_colors)
    rowval = system.generic_matrix.rowval
    colptr = system.generic_matrix.colptr
    nzval = system.generic_matrix.nzval
    for i = 1:(length(colptr) - 1)
        for j = colptr[i]:(colptr[i + 1] - 1)
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

################################################################

mutable struct FactorizationPreconditioner{C}
    cache::C
end

"""
   factorize(A,method::LinearSolve.AbstractFactorization)

Calculate an LU factorization of A using one of the methods available in LinearSolve.jl. 
"""
function LinearAlgebra.factorize(A, method::LinearSolve.AbstractFactorization)
    pr = LinearProblem(A, zeros(eltype(A), size(A, 1)))
    p = FactorizationPreconditioner(init(pr, method))
    p
end

function LinearAlgebra.ldiv!(u, p::FactorizationPreconditioner, b)
    p.cache = LinearSolve.set_b(p.cache, b)
    sol = solve(p.cache)
    p.cache = sol.cache
    copyto!(u, sol.u)
end

function _solve_linear!(u, system, nlhistory, control, method_linear, A, b)
    if isnothing(system.linear_cache)
        Pl = control.precon_linear(A)
        p = LinearProblem(A, b)
        system.linear_cache = init(p,
                                   method_linear;
                                   abstol = control.abstol_linear,
                                   reltol = control.reltol_linear,
                                   verbose = doprint(control, 'l'),
                                   Pl)
    else
        system.linear_cache = LinearSolve.set_A(system.linear_cache, A)
        system.linear_cache = LinearSolve.set_b(system.linear_cache, b)
        if control.keepcurrent_linear
            Pl = control.precon_linear(A)
            system.linear_cache = LinearSolve.set_prec(system.linear_cache, Pl, LinearSolve.Identity())
        end
    end

    try
        sol = LinearSolve.solve(system.linear_cache)
        system.linear_cache = sol.cache
        u .= sol.u
        nliniter = sol.iters
        nlhistory.nlu += 1
        nlhistory.nlin = sol.iters
    catch err
        if (control.handle_exceptions)
            _print_error(err, stacktrace(catch_backtrace()))
            throw(LinearSolverError())
        else
            rethrow(err)
        end
    end
end

"""
$(SIGNATURES)

Solve time step problem. This is the core routine
for implicit Euler and stationary solve
"""
function _solve_timestep!(solution::AbstractMatrix{Tv}, # old time step solution resp. initial value
                          oldsol::AbstractMatrix{Tv}, # old time step solution resp. initial value
                          system::AbstractSystem{Tv, Tc, Ti, Tm}, # Finite volume system
                          control::SolverControl,
                          time,
                          tstep,
                          embedparam,
                          params;
                          called_from_API = false) where {Tv, Tc, Ti, Tm}
    _complete!(system; create_newtonvectors = true)
    nlhistory = NewtonSolverHistory()
    t = @elapsed begin
        solution .= oldsol
        residual = system.residual
        update = system.update
        _initialize!(solution, system; time, λ = embedparam, params)

        method_linear = control.method_linear
        if isnothing(method_linear)
            if Tv != Float64
                method_linear = SparspakFactorization()
            else
                if dim_grid(system.grid) == 1
                    method_linear = KLUFactorization()
                elseif dim_grid(system.grid) == 2
                    method_linear = SparspakFactorization()
                else
                    method_linear = UMFPACKFactorization()
                end
            end
        end

        oldnorm = 1.0
        converged = false
        damp = 1.0
        if !system.is_linear
            if doprint(control, 'n')
                println("  [n]ewton: #it(lin) |update|  cont3tion |round|  #round\n")
            end
            damp = control.damp_initial
            rnorm = control.rnorm(solution)
        end

        nlu_reuse = 0
        nround = 0
        tolx = 0.0
        ncalloc = 0
        nballoc = 0
        neval = 0
        niter = 1

        while niter <= control.maxiters
            # Create Jacobi matrix and RHS for Newton iteration
            try
                nca, nba, nev = eval_and_assemble(system,
                                                  solution,
                                                  oldsol,
                                                  residual,
                                                  time,
                                                  tstep,
                                                  embedparam,
                                                  params;
                                                  edge_cutoff = control.edge_cutoff)
                ncalloc += nca
                nballoc += nba
                neval += nev
            catch err
                if (control.handle_exceptions)
                    _print_error(err, stacktrace(catch_backtrace()))
                    throw(AssemblyError())
                else
                    rethrow(err)
                end
            end

            _solve_linear!(values(update),
                           system,
                           nlhistory,
                           control,
                           method_linear,
                           system.matrix,
                           values(residual))

            values(solution) .-= damp * values(update)

            if system.is_linear
                converged = true
                break
            end

            damp = min(damp * control.damp_growth, 1.0)
            norm = control.unorm(update)
            if tolx == 0.0
                tolx = norm * control.reltol
            end
            dnorm = 1.0
            rnorm_new = control.rnorm(solution)
            if rnorm > 1.0e-50
                dnorm = abs((rnorm - rnorm_new) / rnorm)
            end

            if dnorm < control.tol_round
                nround = nround + 1
            else
                nround = 0
            end

            if control.log
                push!(nlhistory.l1normdiff, dnorm)
                push!(nlhistory.updatenorm, norm)
            end
            if doprint(control, 'n')
                if control.reltol_linear < 1.0
                    itstring = @sprintf("  [n]ewton: % 3d(% 3d)", niter, nlhistory.nlin)
                else
                    itstring = @sprintf("it=% 3d", niter)
                end
                if control.max_round > 0
                    @printf("%s %.3e %.3e %.3e % 2d\n",
                            itstring,
                            norm,
                            norm/oldnorm,
                            dnorm,
                            nround)
                else
                    @printf("%s %.3e %.3e\n", itstring, norm, norm/oldnorm)
                end
            end
            if niter > 1 && norm / oldnorm > 1.0 / control.tol_mono
                converged = false
                break
            end

            if norm < control.abstol || norm < tolx
                converged = true
                break
            end
            oldnorm = norm
            rnorm = rnorm_new

            if nround > control.max_round
                converged = true
                break
            end
            niter = niter + 1
        end
        if !converged
            throw(ConvergenceError())
        end
    end
    if control.log
        nlhistory.time = t
    end

    if ncalloc + nballoc > 0 && doprint(control, 'a')
        @warn "[a]llocations in assembly loop: cells: $(ncalloc÷neval), bfaces: $(nballoc÷neval)"
    end

    if doprint(control, 'n') && !system.is_linear
        println("  [n]ewton: $(round(t,sigdigits=3)) seconds")
    end

    if doprint(control, 'l') && system.is_linear
        println("  [l]inear($(nameof(typeof(method_linear)))): $(round(t,sigdigits=3)) seconds")
    end

    system.history = nlhistory
end

################################################################
"""
````
solve!(solution, inival, system; 
    control=SolverControl(), 
    tstep=Inf)
````
Mutating version of [`solve(inival,system)`](@ref)
"""
function VoronoiFVM.solve!(solution, # Solution
                           inival,   # Initial value 
                           system::VoronoiFVM.AbstractSystem;     # Finite volume system
                           control = SolverControl(),      # Newton solver control information
                           time = Inf,
                           tstep = Inf,                # Time step size. Inf means  stationary solution
                           embedparam = 0.0,
                           params = zeros(0),
                           called_from_API = false)
    fix_deprecations!(control)
    if !called_from_API && doprint(control, 'd')
        @warn "[d]eprecated: solve(inival,solution, system; kwargs...)"
    end
    _solve_timestep!(solution,
                     inival,
                     system,
                     control,
                     time,
                     tstep,
                     embedparam,
                     params;
                     called_from_API = true)
    return solution
end

################################################################
"""
````
    solve(inival, system; control=SolverControl(),params, tstep=Inf)
````
Alias for [`solve(system::VoronoiFVM.AbstractSystem; kwargs...)`](@ref) with the corresponding keyword arguments.

Solve stationary problem(if `tstep==Inf`) or one step implicit Euler step using Newton's method with `inival` as initial
value. Returns a solution array.
"""
function CommonSolve.solve(inival,   # Initial value 
                           system::AbstractSystem;     # Finite volume system
                           control = SolverControl(),      # Newton solver control information
                           time = Inf,
                           tstep = Inf,                # Time step size. Inf means  stationary solution
                           params = zeros(0),
                           called_from_API = false)
    fix_deprecations!(control)
    if !called_from_API && doprint(control, 'd')
        @warn "[d]eprecated: solve(inival,system; kwargs...)"
    end

    solve!(unknowns(system),
           inival,
           system;
           control = control,
           time = time,
           tstep = tstep,
           params = params,
           called_from_API = true)
end

Δλ_val(control, transient) = transient ? control.Δt : control.Δp
Δλ_min(control, transient) = transient ? control.Δt_min : control.Δp_min
Δλ_max(control, transient) = transient ? control.Δt_max : control.Δp_max
Δλ_grow(control, transient) = transient ? control.Δt_grow : control.Δp_grow

"""
        solve(inival, system, times; kwargs...)

Alias for [`solve(system::VoronoiFVM.AbstractSystem; kwargs...)`](@ref) with the corresponding keyword arguments.

"""
function CommonSolve.solve(inival,
                           system::VoronoiFVM.AbstractSystem,
                           lambdas;
                           control = SolverControl(),
                           transient = true, # choose between transient and stationary (embedding) case
                           time = 0.0,
                           params = zeros(0),
                           called_from_API = false,
                           kwargs...)
    fix_deprecations!(control)
    if !called_from_API && doprint(control, 'd')
        @warn "[d]eprecated: solve(inival,system,times;kwargs...)"
    end
    λstr = "t"
    if !transient
        λstr = "p"
    end

    allhistory = TransientSolverHistory()

    solution = copy(inival)
    oldsolution = copy(inival)
    _initialize_dirichlet!(inival, system; time, λ = Float64(lambdas[1]), params)
    Δλ = Δλ_val(control, transient)

    t0 = @elapsed if !transient
        control.pre(solution, Float64(lambdas[1]))
        solution = solve!(solution,
                          oldsolution,
                          system;
                          called_from_API = true,
                          control = control,
                          time = time,
                          tstep = Inf,
                          embedparam = Float64(lambdas[1]),
                          params = params,
                          kwargs...)
        control.post(solution, oldsolution, lambdas[1], 0)
        if control.log
            push!(allhistory, system.history)
            push!(allhistory.times, lambdas[1])
            Δu = control.delta(system, solution, oldsolution, lambdas[1], 0)
            push!(allhistory.updates, Δu)
        end
        oldsolution .= solution
    end

    tsol = TransientSolution(Float64(lambdas[1]), solution; in_memory = control.in_memory)

    if doprint(control, 'e')
        println("[e]volution: start in $(extrema(lambdas))")
    end

    istep = 0
    t1 = @elapsed for i = 1:(length(lambdas) - 1)
        Δλ = max(Δλ, Δλ_min(control, transient))
        λstart = lambdas[i]
        λend = lambdas[i + 1]
        λ = Float64(λstart)

        while λ < λend
            solved = false
            λ0 = λ
            Δu = 0.0
            while !solved
                solved = true
                try
                    λ = λ0 + Δλ
                    control.pre(solution, λ)
                    if transient
                        solution = solve!(solution,
                                          oldsolution,
                                          system;
                                          called_from_API = true,
                                          control = control,
                                          time = λ,
                                          tstep = Δλ,
                                          params = params,
                                          kwargs...)
                    else
                        solution = solve!(solution,
                                          oldsolution,
                                          system;
                                          called_from_API = true,
                                          control = control,
                                          time = time,
                                          tstep = Inf,
                                          embedparam = λ,
                                          params = params,
                                          kwargs...)
                    end
                catch err
                    if (control.handle_exceptions)
                        _print_error(err, stacktrace(catch_backtrace()))
                    else
                        rethrow(err)
                    end
                    solved = false
                end
                if solved
                    Δu = control.delta(system, solution, oldsolution, λ, Δλ)
                    if Δu > 2.0 * control.Δu_opt
                        solved = false
                    end
                end
                if !solved
                    # reduce time step and retry  solution
                    Δλ = Δλ * 0.5
                    if Δλ < Δλ_min(control, transient)
                        if !(control.force_first_step && istep == 0)
                            throw(EmbeddingError(" Δ$(λstr)_min=$(Δλ_min(control,transient)) reached while Δu=$(Δu) >>  Δu_opt=$(control.Δu_opt) "))
                        else
                            solved = true
                        end
                    end
                    if doprint(control, 'e')
                        @printf("[e]volution:  Δu=%.3e => retry: Δ%s=%.3e\n", Δu, λstr, Δλ)
                    end
                end
            end
            istep = istep + 1
            if doprint(control, 'e')
                @printf("[e]volution: step=%d %s=%.3e Δ%s=%.3e Δu=%.3e\n",
                        istep,
                        λstr,
                        λ,
                        λstr,
                        Δλ,
                        Δu)
            end
            if control.log
                push!(allhistory, system.history)
                push!(allhistory.updates, Δu)
                push!(allhistory.times, λ)
            end
            if control.store_all
                append!(tsol, λ, solution)
            end
            control.post(solution, oldsolution, λ, Δλ)
            oldsolution .= solution
            if λ < λend
                Δλ = min(Δλ_max(control, transient),
                         Δλ * Δλ_grow(control, transient),
                         Δλ * control.Δu_opt / (Δu + 1.0e-14),
                         λend - λ)
            end
        end
        if !control.store_all
            append!(tsol, lambdas[i + 1], solution)
        end
        control.sample(solution, lambdas[i + 1])
    end

    if doprint(control, 'e')
        println("[e]volution:  $(round(t0+t1,sigdigits=3)) seconds")
    end

    system.history = allhistory
    return tsol
end

"""
    evaluate_residual_and_jacobian(system,u;
                                   t=0.0, tstep=Inf,embed=0.0)

Evaluate residual and jacobian at solution value u.
Returns a solution vector containing the residual, and an ExendableSparseMatrix
containing the linearization at u.

"""
function evaluate_residual_and_jacobian(sys, u; t = 0.0, tstep = Inf, embed = 0.0)
    _complete!(sys; create_newtonvectors = true)

    eval_and_assemble(sys, u, u, sys.residual, t, tstep, embed, zeros(0))
    copy(sys.residual), copy(flush!(sys.matrix))
end

module NoModule end

#####################################################################
"""
    solve(system; kwargs...)
    
Solution method for [`VoronoiFVM.System`](@ref).  
For using ODE solvers from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), see
the [VoronoiFVMDiffEq.jl](https://github.com/j-fu/VoronoiFVMDiffEq.jl) package.
    
Keyword arguments:
- General for all solvers 
   - `inival` (default: 0) : Array created via [`unknowns`](@ref) or  number giving the initial value.
   -  All elements of [`SolverControl`](@ref) can be used as kwargs.
   - `control` (default: nothing): Pass instance of [`SolverControl`](@ref)
   - `params`: Parameters (Parameter handling is experimental and may change)
    
- __Stationary solver__:
  Invoked if neither `times` nor `embed`, nor `tstep` are given as keyword argument.
  - `time` (default: `0.0`): Set time value.
  Returns a [`DenseSolutionArray`](@ref) or [`SparseSolutionArray`](@ref)

- __Embedding (homotopy) solver__: Invoked if `embed` kwarg is given.
  Use homotopy embedding + damped Newton's method  to 
  solve stationary problem or to solve series of parameter dependent problems.
  Parameter step control is performed according to solver control data.  kwargs and default values are:
  - `embed` (default: `nothing` ): vector of parameter values to be reached exactly
  In addition,  all kwargs of the implicit Euler solver (besides `times`) are handled.  
  Returns a transient solution object `sol` containing the stored solution(s),  see [`TransientSolution`](@ref).
  
- __Implicit Euler transient solver__: Invoked if `times` kwarg is given.
  Use implicit Euler method  + damped   Newton's method  to 
  solve time dependent problem. Time step control is performed
  according to solver control data.  kwargs and default values are:
  - `times` (default: `nothing` ): vector of time values to be reached exactly
  - `pre` (default: `(sol,t)->nothing` ):  invoked before each time step
  - `post`  (default:  `(sol,oldsol, t, Δt)->nothing` ):  invoked after each time step
  - `sample` (default:  `(sol,t)->nothing` ): invoked after timestep for all times in `times[2:end]`.
  - `delta` (default:  `(system, u,v,t, Δt)->norm(sys,u-v,Inf)` ):  Value  used to control the time step size `Δu`
  If `control.handle_error` is true, if step solution  throws an error,
  stepsize  is lowered, and  step solution is called again with a smaller time value.
  If `control.Δt<control.Δt_min`, solution is aborted with error.
  Returns a transient solution object `sol` containing the stored solution,  see [`TransientSolution`](@ref).
  
- __Implicit Euler timestep solver__.  Invoked if `tstep` kwarg is given. Solve one time step of the implicit Euler method.
  - `time` (default: `0`): Set time value. 
  - `tstep`: time step
  Returns a [`DenseSolutionArray`](@ref) or [`SparseSolutionArray`](@ref)
"""
function CommonSolve.solve(sys::VoronoiFVM.AbstractSystem;
                           inival = 0,
                           params = zeros(0),
                           control = VoronoiFVM.SolverControl(),
                           time = 0.0,
                           tstep = Inf,
                           kwargs...)
    fix_deprecations!(control)

    if isa(inival, Number)
        inival = unknowns(sys; inival = inival)
    elseif !VoronoiFVM.isunknownsof(inival, sys)
        @error "wrong type of inival: $(typeof(inival))"
    end

    for pair in kwargs
        if first(pair) != :times &&
           first(pair) != :embed &&
           hasfield(SolverControl, first(pair))
            setproperty!(control, first(pair), last(pair))
        end
    end

    sys.linear_cache = nothing

    if haskey(kwargs, :times) && !isnothing(kwargs[:times])
        solve(inival,
              sys,
              kwargs[:times];
              control,
              transient = true,
              params,
              time = kwargs[:times][1],
              called_from_API = true)
    elseif haskey(kwargs, :embed) && !isnothing(kwargs[:embed])
        solve(inival,
              sys,
              kwargs[:embed];
              called_from_API = true,
              transient = false,
              control,
              params,
              time)
    else
        solve(inival, sys; called_from_API = true, control, params, time, tstep)
    end
end
