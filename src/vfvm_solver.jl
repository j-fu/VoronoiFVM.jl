##################################################################
"""
$(SIGNATURES)

Extract value from dual number. Use to debug physics callbacks.
Re-exported from ForwardDiff.jl
"""
const value=ForwardDiff.value

# These are needed to enable iterative solvers to work with dual numbers
Base.Float64(x::ForwardDiff.Dual)=value(x)
Random.rand(rng::AbstractRNG, ::Random.SamplerType{ForwardDiff.Dual{T,V,N}}) where {T,V,N} = ForwardDiff.Dual{T,V,N}(rand(rng,V))


##################################################################
"""
$(SIGNATURES)

Extract derivatives from dual number.
"""
const partials=ForwardDiff.partials

##################################################################
"""
$(SIGNATURES)

Extract number of derivatives from dual number.
"""
const npartials=ForwardDiff.npartials


"""
Add value to matrix if it is nonzero
"""
@inline function _addnz(matrix,i,j,v::Tv,fac) where Tv
    if isnan(v)
        error("trying to assemble NaN")
    end
    if v!=zero(Tv)
        rawupdateindex!(matrix,+,v*fac,i,j)
    end
end

ExtendableSparse.rawupdateindex!(m::AbstractMatrix,op,v,i,j)= m[i,j]=op(m[i,j],v)


"""
$(TYPEDEF)

Exception thrown if Newton's method convergence fails.
"""
struct ConvergenceError <: Exception
end

"""
$(TYPEDEF)

Exception thrown if error occured during assembly (e.g. domain error)
"""
struct AssemblyError <: Exception
end

"""
$(TYPEDEF)

Exception thrown if error occured during factorization.
"""
struct FactorizationError <: Exception
end

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
function _print_error(err,st)
    println()
    println(err)
    nlines=5
    for i=1:min(nlines,length(st))
        line=@sprintf("%s",st[i])
        L=length(line)
        if L<80
            println(line)
        else
            print(line[1:35])
            print(" ... ")
            println(line[L-35:L])
        end
    end
    if length(st)>nlines
        println("...")
    end
    println()
end


"""
Solve time step problem. This is the core routine
for implicit Euler and stationary solve
"""
function _solve!(
    solution::AbstractMatrix{Tv}, # old time step solution resp. initial value
    oldsol::AbstractMatrix{Tv}, # old time step solution resp. initial value
    system::AbstractSystem{Tv, Tc, Ti, Tm}, # Finite volume system
    control::NewtonControl,
    time,
    tstep,
    embedparam,
    params;
    mynorm=(u)->LinearAlgebra.norm(values(u),Inf),
    myrnorm=(u)->LinearAlgebra.norm(values(u),1),
    kwargs...
) where {Tv, Tc, Ti, Tm}

    _complete!(system, create_newtonvectors=true)
    nlhistory=NewtonSolverHistory()
    t=@elapsed begin
        solution.=oldsol
        residual=system.residual
        update=system.update
        _initialize!(solution,system; time, λ=embedparam,params)
        if VERSION<=v"1.8"
            SuiteSparse.UMFPACK.umf_ctrl[3+1]=control.umfpack_pivot_tolerance
        end
        
        if isnothing(system.factorization)
            if isa(system.matrix,ExtendableSparseMatrix)
                system.factorization=factorization(control; valuetype=Tv)
            end
        end
        
        oldnorm=1.0
        converged=false
        if control.verbose
            @printf("    Start Newton iteration\n")
        end
        nlu_reuse=0
        nround=0
        damp=control.damp_initial
        tolx=0.0
        rnorm=myrnorm(solution)

        for ii=1:control.max_iterations
            # Create Jacobi matrix and RHS for Newton iteration
            try
                eval_and_assemble(system,solution,oldsol,residual,time,tstep,embedparam,params,edge_cutoff=control.edge_cutoff)
            catch err
                if (control.handle_exceptions)
                    _print_error(err,stacktrace(catch_backtrace()))
                    throw(AssemblyError())
                else
                    rethrow(err)
                end
            end
            
            
            ## Linear system solution: mtx*update=residual
            mtx=system.matrix
            nliniter=0
            
            ## (re)factorize, if possible reusing old factorization data
            if !isnothing(system.factorization) && nlu_reuse==0 
                try
                    factorize!(system.factorization,mtx)
                catch err
                    if (control.handle_exceptions)
                        _print_error(err,stacktrace(catch_backtrace()))
                        throw(FactorizationError())
                    else
                        rethrow(err)
                    end
                end
                if control.log
                    nlhistory.nlu+=1
                end
            end
            
            ## Direct solution via LU solve:
            if !isnothing(system.factorization)  && issolver(system.factorization) && nlu_reuse==0
                try
                    ldiv!(values(update),system.factorization,values(residual))
                    if control.log
                        nlhistory.nlin+=1
                    end
                catch err
                    if (control.handle_exceptions)
                        _print_error(err,stacktrace(catch_backtrace()))
                        throw(FactorizationError())
                    else
                        rethrow(err)
                    end
                end
                
            ## Iterative solution using factorization as preconditioner
            elseif !isnothing(system.factorization) && !issolver(system.factorization)  || nlu_reuse>0
                if control.iteration == :krylov_bicgstab
                    (sol,history)= Krylov.bicgstab(mtx,
                                                   values(residual),
                                                   M=system.factorization,
                                                   rtol=control.tol_linear,
                                                   ldiv=true,
                                                   itmax=control.max_iterations_linear,
                                                   history=true)
                    values(update).=sol
                elseif control.iteration == :krylov_gmres
                    (sol,history)= Krylov.gmres(mtx,
                                                values(residual),
                                                M=system.factorization,
                                                rtol=control.tol_linear,
                                                ldiv=true,
                                                restart=true,
                                                itmax=control.max_iterations_linear,
                                                memory=control.gmres_restart,
                                                history=true)
                    values(update).=sol
                elseif control.iteration == :krylov_cg
                    (sol,history)= Krylov.cg(mtx,
                                             values(residual),
                                             M=system.factorization,
                                             ldiv=true,
                                             rtol=control.tol_linear,
                                             itmax=control.max_iterations_linear,
                                             history=true)
                    values(update).=sol
                elseif control.iteration == :cg
                    update.=zero(Tv)
                    (sol,history)= IterativeSolvers.cg!(values(update),
                                                        mtx,
                                                        values(residual),
                                                        Pl=system.factorization,
                                                        reltol=control.tol_linear,
                                                        maxiter=control.max_iterations_linear,
                                                        log=true)
                elseif control.iteration == :bicgstab
                    update.=zero(Tv)
                    (sol,history)= IterativeSolvers.bicgstabl!(values(update),
                                                               mtx,
                                                               values(residual),
                                                               1,
                                                               Pl=system.factorization,
                                                               reltol=control.tol_linear,
                                                               max_mv_products=control.max_iterations_linear,
                                                               log=true)
                else
                    error("wrong value of `iteration`, see documentation of SolverControl")
                end

                if control.iteration == :cg || control.iteration == :bicgstab 
                    nliniter=history.iters
                else
                    nliniter=history.niter
                end
                if control.log
                    nlhistory.nlin+=history.niter
                end
            else
                values(update).=system.matrix\values(residual)
            end

            if control.max_lureuse>0
                nlu_reuse=(nlu_reuse+1)%control.max_lureuse
            end
            solval=values(solution)
            solval.-=damp*values(update)
            damp=min(damp*control.damp_growth,1.0)
            norm=mynorm(update)
            if tolx==0.0
                tolx=norm*control.tol_relative
            end
            dnorm=1.0
            rnorm_new=myrnorm(solution)
            if rnorm>1.0e-50
                dnorm=abs((rnorm-rnorm_new)/rnorm)
            end
            
            if dnorm<control.tol_round
                nround=nround+1
            else
                nround=0
            end

            if control.log
                push!(nlhistory.l1normdiff,dnorm)
                push!(nlhistory.updatenorm,norm)
            end
            if control.verbose
                if   control.tol_linear<1.0
                itstring=@sprintf("it=% 3d(% 2d)",ii,nliniter)
                else
                    itstring=@sprintf("it=% 3d",ii)
            end
                if control.max_round>0
                    @printf("    %s du=%.3e cont=%.3e dnorm=%.3e %d\n",itstring,norm, norm/oldnorm,dnorm,nround)
                else
                    @printf("    %s du=%.3e cont=%.3e\n",itstring,norm, norm/oldnorm)
            end
            end
            if ii>1 &&  norm/oldnorm > 1.0/control.tol_mono
                converged=false
            break
            end
            
            if norm<control.tol_absolute || norm <tolx
                converged=true
                break
            end
            oldnorm=norm
            rnorm=rnorm_new
            
            if nround>control.max_round
            converged=true
                break
            end
        end
        if !converged
            throw(ConvergenceError())
        end
    end
    if control.log
        nlhistory.time=t
    end
    if control.verbose
       @printf("    Newton iteration successful\n")
    end
    system.history=nlhistory
end

function zero!(m::ExtendableSparseMatrix{Tv,Ti}) where {Tv,Ti}
    nzv  = nonzeros(m)
    nzv .= zero(Tv)
end

zero!(m::AbstractMatrix{T}) where T = m.=zero(T)

# function zero!(m::MultidiagonalMatrix{T}) where {T}
#     for d ∈ m.shadow
#         d.second.=0.0
#     end
# end

################################################################
"""
$(SIGNATURES)

Main assembly method.

Evaluate solution with result in right hand side F and 
assemble Jacobi matrix into system.matrix.
"""
function eval_and_assemble(system::AbstractSystem{Tv, Tc, Ti, Tm},
                           U::AbstractMatrix{Tv}, # Actual solution iteration
                           UOld::AbstractMatrix{Tv}, # Old timestep solution
                           F::AbstractMatrix{Tv},# Right hand side
                           time,
                           tstep,# time step size. Inf means stationary solution
                           λ,
                           params::AbstractVector;
                           edge_cutoff=0.0
                           ) where {Tv, Tc, Ti, Tm}
    

    _complete!(system) # needed here as well for test function system which does not use newton
    
    grid    = system.grid
    physics = system.physics
    node    = Node(system,time,λ,params)
    bnode   = BNode(system,time,λ,params)
    edge    = Edge(system,time,λ,params)
    bedge   = BEdge(system,time,λ,params)

    
    
    nspecies::Int = num_species(system)
    matrix        = system.matrix
    
    cellnodefactors::Array{Tv,2}  = system.cellnodefactors
    celledgefactors::Array{Tv,2}  = system.celledgefactors
    bfacenodefactors::Array{Tv,2} = system.bfacenodefactors
    bfaceedgefactors::Array{Tv,2} = system.bfaceedgefactors    
    
    # Reset matrix + rhs
    zero!(matrix)
    F   .= 0.0
    nparams::Int=system.num_parameters
    if nparams>0
        dudp=system.dudp
    end
    for iparam=1:nparams
        dudp[iparam].=0.0
    end

    # Arrays for gathering solution data
    UK    = Array{Tv,1}(undef,nspecies + nparams)
    UKOld = Array{Tv,1}(undef,nspecies + nparams)
    UKL   = Array{Tv,1}(undef,2*nspecies + nparams)

    @assert length(params)==nparams
    if nparams>0
        UK[nspecies+1:end].=params
        UKOld[nspecies+1:end].=params
        UKL[2*nspecies+1:end].=params
    end
    
    
    # Inverse of timestep
    # According to Julia documentation, 1/Inf=0 which
    # comes handy to write compact code here for the
    # case of stationary problems.
    tstepinv = 1.0/tstep 

    
    boundary_factors::Array{Tv,2} = system.boundary_factors
    boundary_values::Array{Tv,2}  = system.boundary_values
    bfaceregions::Vector{Ti} = grid[BFaceRegions]
    cellregions::Vector{Ti} = grid[CellRegions]

    nbfaces = num_bfaces(grid)
    ncells  = num_cells(grid)
    geom=grid[CellGeometries][1]
    bgeom   = grid[BFaceGeometries][1]
    
    nn::Int = num_nodes(geom)
    ne::Int = num_edges(geom)

    has_legacy_bc = !iszero(boundary_factors) || !iszero(boundary_values)


    #
    # These wrap the different physics functions.
    #
    src_eval  = ResEvaluator(physics,:source,UK,node,nspecies)
    rea_eval  = ResJacEvaluator(physics,:reaction,UK,node,nspecies)
    stor_eval = ResJacEvaluator(physics,:storage,UK,node,nspecies)
    oldstor_eval = ResEvaluator(physics,:storage,UK,node,nspecies)
    flux_eval = ResJacEvaluator(physics,:flux,UKL,edge,nspecies)


    bsrc_eval  = ResEvaluator(physics,:bsource,UK,bnode,nspecies)
    brea_eval  = ResJacEvaluator(physics,:breaction,UK,bnode,nspecies)
    bstor_eval = ResJacEvaluator(physics,:bstorage,UK,bnode,nspecies)
    oldbstor_eval = ResEvaluator(physics,:bstorage,UK,bnode,nspecies)
    bflux_eval = ResJacEvaluator(physics,:bflux,UKL,bedge,nspecies)


    ncalloc=@allocated  for icell=1:ncells
        ireg=cellregions[icell]
        for inode=1:nn
            _fill!(node,inode,icell)
            @views UK[1:nspecies]    .= U[:,node.index]
            @views UKOld[1:nspecies] .= UOld[:,node.index]
            # xx gather:
                # ii=0
                # for i region_spec.colptr[ireg]:region_spec.colptr[ireg+1]-1
                #    ispec=Fdof.rowval[idof]
                #    ii=ii+1
                #    node.speclist[ii]=ispec
                #    UK[ii]=U[ispec,K]
                #   UK[1:nspecies]=U[:,node.index]
                #   UKOld[1:nspecies]=UOld[:,node.index]
                # Evaluate source term

            evaluate!(src_eval)
            src = res(src_eval)
            
            evaluate!(rea_eval,UK)
            res_react = res(rea_eval)
            jac_react = jac(rea_eval)
            
            evaluate!(stor_eval,UK)
            res_stor = res(stor_eval)
            jac_stor = jac(stor_eval)

            evaluate!(oldstor_eval,UKOld)
            oldstor = res(oldstor_eval)

            K=node.index
            for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                ispec=_spec(F,idof,K)
                if system.region_species[ispec,ireg]<=0
                    continue
                end
                _add(F,idof,cellnodefactors[inode,icell]*(res_react[ispec]-src[ispec] + (res_stor[ispec]-oldstor[ispec])*tstepinv))
                for jdof=_firstnodedof(F,K):_lastnodedof(F,K)
                    jspec=_spec(F,jdof,K)
                    if system.region_species[jspec,ireg]<=0
                        continue
                    end
                    val=jac_react[ispec,jspec]+ jac_stor[ispec,jspec]*tstepinv
                    _addnz(matrix,idof,jdof,val,cellnodefactors[inode,icell])
                end
                for iparam=1:nparams
                    jparam=nspecies+iparam
                    dudp[iparam][idof]+=(jac_react[ispec,jparam]+ jac_stor[ispec,jparam]*tstepinv)*cellnodefactors[inode,icell]
                end
            end
        end

        for iedge=1:ne
            if abs(celledgefactors[iedge,icell])<edge_cutoff
                continue
            end
            _fill!(edge,iedge,icell)

            #Set up argument for fluxwrap
            @views UKL[1:nspecies]            .= U[:, edge.node[1]]
            @views UKL[nspecies+1:2*nspecies] .= U[:, edge.node[2]]


            evaluate!(flux_eval,UKL)
            res_flux = res(flux_eval)
            jac_flux = jac(flux_eval)

            K=edge.node[1]
            L=edge.node[2]
            
            fac=celledgefactors[iedge,icell]
            for idofK=_firstnodedof(F,K):_lastnodedof(F,K)
                ispec=_spec(F,idofK,K)
                idofL=dof(F,ispec,L)
                if idofL==0
                    continue
                end
                if system.region_species[ispec,ireg]<=0
                     continue
                end

                _add(F,idofK,fac*res_flux[ispec])
                _add(F,idofL,-fac*res_flux[ispec])
                
                for jdofK=_firstnodedof(F,K):_lastnodedof(F,K)
                    jspec=_spec(F,jdofK,K)
                    if system.region_species[jspec,ireg]<=0
                        continue
                    end
                    jdofL=dof(F,jspec,L)
                    if jdofL==0
                        continue
                    end
                    
                    _addnz(matrix,idofK,jdofK,+jac_flux[ispec,jspec         ],fac)
                    _addnz(matrix,idofL,jdofK,-jac_flux[ispec,jspec         ],fac)
                    _addnz(matrix,idofK,jdofL,+jac_flux[ispec,jspec+nspecies],fac)
                    _addnz(matrix,idofL,jdofL,-jac_flux[ispec,jspec+nspecies],fac)
                end
                
                for iparam=1:nparams
                    jparam=2*nspecies+iparam
                    dudp[iparam][idofK]+=fac*jac_flux[ispec,jparam]
                    dudp[iparam][idofL]-=fac*jac_flux[ispec,jparam]
                end
            end
        end
    end

    nbn::Int = num_nodes(bgeom)
    nbe::Int = num_edges(bgeom)

    # Assembly loop for boundary conditions
    nballoc= @allocated for ibface=1:nbfaces
        ibreg=bfaceregions[ibface]

        # Loop over nodes of boundary face
        for ibnode=1:nbn
            # Fill bnode data shuttle with data from grid
            _fill!(bnode,ibnode,ibface)

            
            # Copy unknown values from solution into dense array
            @views UK[1:nspecies].=U[:,bnode.index]

            # Measure of boundary face part assembled to node
            bnode_factor::Tv=bfacenodefactors[ibnode,ibface]
            bnode.Dirichlet=Dirichlet/bnode_factor
            
            evaluate!(bsrc_eval)
            bsrc = res(bsrc_eval)

            # if isbsource
            #     bsourcewrap(bsrc)
            # end

            # Global index of node
            K::Int=bnode.index


            if has_legacy_bc
                # Assemble "standard" boundary conditions: Robin or
                # Dirichlet
                # valid only for interior species, currently not checked
                for ispec=1:nspecies
                    idof=dof(F,ispec,K)
                    # If species is present, assemble the boundary condition
                    if idof>0
                        # Get user specified data
                        boundary_factor = boundary_factors[ispec,ibreg]
                        boundary_value  = boundary_values[ispec,ibreg]
                        
                        if boundary_factor==Dirichlet
                            # Dirichlet is encoded in the boundary condition factor
                            # Penalty method: scale   (u-boundary_value) by penalty=boundary_factor
                            
                            # Add penalty*boundry_value to right hand side
                            F[ispec,K]+=boundary_factor*(U[ispec,K]-boundary_value)
                            
                            # Add penalty to matrix main diagonal (without bnode factor, so penalty
                            # is independent of h)
                            _addnz(matrix,idof,idof,boundary_factor,1)
                        else
                            # Robin boundary condition
                            F[ispec,K]+=bnode_factor*(boundary_factor*U[ispec,K]-boundary_value)
                            _addnz(matrix,idof,idof,boundary_factor, bnode_factor)
                        end
                    end
                end
            end

            if docall(brea_eval)
                evaluate!(brea_eval,UK)
                res_breact = res(brea_eval)
                jac_breact = jac(brea_eval)
                
                # Assemble RHS + matrix
                for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                    ispec=_spec(F,idof,K)
                    if !isdof(system,ispec,K)
                        continue
                    end
                    
                    _add(F,idof,bnode_factor*(res_breact[ispec]-bsrc[ispec]))
                    for jdof=_firstnodedof(F,K):_lastnodedof(F,K)
                        jspec=_spec(F,jdof,K)
                        if !isdof(system,jspec,K)
                            continue
                        end
                        _addnz(matrix,idof,jdof,jac_breact[ispec,jspec],bnode_factor)
                    end
                    for iparam=1:nparams
                        jparam=nspecies+iparam
                        dudp[iparam][idof]+=jac_breact[ispec,jparam]*bnode_factor
                    end
                end
            end

            if docall(bstor_eval)
                evaluate!(bstor_eval,UK)
                res_bstor = res(bstor_eval)
                jac_bstor = jac(bstor_eval)
                
                @views UKOld.=UOld[:,bnode.index]
                evaluate!(oldbstor_eval,UKOld)
                oldbstor = res(oldbstor_eval)
                
                for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                    ispec=_spec(F,idof,K)
                    if !isdof(system,ispec,K)
                        continue
                    end
                    
                    # Assemble finite difference in time for right hand side
                    _add(F,idof,bnode_factor*(res_bstor[ispec]-oldbstor[ispec])*tstepinv)
                    # Assemble matrix.
                    for jdof=_firstnodedof(F,K):_lastnodedof(F,K)
                        jspec=_spec(F,jdof,K)
                        if !isdof(system,jspec,K)
                            continue
                        end
                        
                        _addnz(matrix,idof,jdof,jac_bstor[ispec,jspec],bnode_factor*tstepinv)
                    end
                    for iparam=1:nparams
                        jparam=nspecies+iparam
                    dudp[iparam][idof]+=jac_bstor[ispec,jparam]*bnode_factor*tstepinv
                    end
                end
            end
        end # ibnode=1:nbn
        
        if docall(bflux_eval)
            for ibedge=1:nbe
                if abs(bfaceedgefactors[ibedge,ibface]) < edge_cutoff
                    continue
                end

                _fill!(bedge,ibedge,ibface)
                @views UKL[1:nspecies]            .= U[:, bedge.node[1]]
                @views UKL[nspecies+1:2*nspecies] .= U[:, bedge.node[2]]

                evaluate!(bflux_eval,UKL)
                res_bflux = res(bflux_eval)
                jac_bflux = jac(bflux_eval)

                K   = bedge.node[1]
                L   = bedge.node[2]
                
                fac = bfaceedgefactors[ibedge, ibface]
                
                for idofK = _firstnodedof(F, K):_lastnodedof(F, K)
                    ispec =_spec(F, idofK, K)
                    if !isdof(system,ispec,K)
                        continue
                    end
                    
                    idofL = dof(F, ispec, L)
                    if idofL == 0
                        continue
                    end
                    _add(F, idofK,  fac*res_bflux[ispec])
                    _add(F, idofL, -fac*res_bflux[ispec])
                    
                    for jdofK = _firstnodedof(F,K):_lastnodedof(F,K)
                        jspec = _spec(F,jdofK,K)
                        if !isdof(system,jspec,K)
                            continue
                        end

                        jdofL = dof(F,jspec,L)
                        if jdofL == 0
                            continue
                        end
                        
                        _addnz(matrix, idofK, jdofK, +jac_bflux[ispec, jspec         ], fac)
                        _addnz(matrix, idofL, jdofK, -jac_bflux[ispec, jspec         ], fac)
                        _addnz(matrix, idofK, jdofL, +jac_bflux[ispec, jspec+nspecies], fac)
                        _addnz(matrix, idofL, jdofL, -jac_bflux[ispec, jspec+nspecies], fac)
                     end
                end
                #!!! Param
                
                
            end
        end
    end

    noallocs(m::ExtendableSparseMatrix)=isnothing(m.lnkmatrix)
    noallocs(m::AbstractMatrix)=false
    # if  no new matrix entries have been created, we should see no allocations
    # in the previous two loops
    if noallocs(matrix) && !_check_allocs(system,ncalloc+nballoc)
        error("""Allocations in assembly loop: cells: $(ncalloc), bfaces: $(nballoc)
                            See the documentation of `check_allocs!` for more information""")
    end
   _eval_and_assemble_generic_operator(system,U,F)
   _eval_and_assemble_inactive_species(system,U,UOld,F)
end


"""
Evaluate and assemble jacobian for generic operator part.
"""
function _eval_and_assemble_generic_operator(system::AbstractSystem,U,F)

    if !has_generic_operator(system)
        return
    end
    generic_operator(f,u)=system.physics.generic_operator(f,u,system)
    vecF=values(F)
    vecU=values(U)
    y=similar(vecF)
    generic_operator(y,vecU)
    vecF.+=y
    forwarddiff_color_jacobian!(system.generic_matrix, generic_operator, vecU, colorvec = system.generic_matrix_colors)
    rowval=system.generic_matrix.rowval
    colptr=system.generic_matrix.colptr
    nzval=system.generic_matrix.nzval
    for i=1:length(colptr)-1
        for j=colptr[i]:colptr[i+1]-1
            updateindex!(system.matrix,+,nzval[j],rowval[j],i)
        end
    end
end


function _eval_generic_operator(system::AbstractSystem,U,F)
    if !has_generic_operator(system)
        return
    end
    generic_operator(f,u)=system.physics.generic_operator(f,u,system)
    vecF=values(F)
    vecU=values(U)
    y=similar(vecF)
    generic_operator(y,vecU)
    vecF.+=y
end

################################################################
"""
````
solve!(solution, inival, system; 
    control=NewtonControl(), 
    tstep=Inf)
````
Mutating version of [`solve(inival,system)`](@ref)
"""
function solve!(
    solution, # Solution
    inival,   # Initial value 
    system::VoronoiFVM.AbstractSystem;     # Finite volume system
    control=NewtonControl(),      # Newton solver control information
    time=Inf,
    tstep=Inf,                # Time step size. Inf means  stationary solution
    embedparam=0.0,
    params=zeros(0),
    kwargs...
)
    if control.verbose
        @time begin
            _solve!(solution,inival,system,control,time,tstep,embedparam,params; kwargs...)
        end
    else 
        _solve!(solution,inival,system,control,time,tstep,embedparam,params; kwargs...)
    end
    return solution
end





################################################################
"""
````
solve(inival, system; control=NewtonControl(),params, tstep=Inf)
````
Alias for [`solve(system::VoronoiFVM.AbstractSystem; kwargs...)`](@ref) with the corresponding keyword arguments.

Solve stationary problem(if `tstep==Inf`) or one step implicit Euler step using Newton's method with `inival` as initial
value. Returns a solution array.
"""
function solve(
    inival,   # Initial value 
    system::AbstractSystem;     # Finite volume system
    control=NewtonControl(),      # Newton solver control information
    time=Inf,
    tstep=Inf,                # Time step size. Inf means  stationary solution
    params=zeros(0),
    kwargs...
)
    solve!(unknowns(system),inival,system; control=control,time=time, tstep=tstep, params=params, kwargs...)
end


Δλ_val(control,transient)= transient ? control.Δt : control.Δp
Δλ_min(control,transient)= transient ? control.Δt_min : control.Δp_min
Δλ_max(control,transient)= transient ? control.Δt_max : control.Δp_max
Δλ_grow(control,transient)= transient ? control.Δt_grow : control.Δp_grow



"""
        solve(inival, system, times; kwargs...)

Alias for [`solve(system::VoronoiFVM.AbstractSystem; kwargs...)`](@ref) with the corresponding keyword arguments.

"""
function solve(inival,
               system::VoronoiFVM.AbstractSystem,
               lambdas;
               control=NewtonControl(),
               pre=function(sol,t) end,       # Function for preparing step
               post=function(sol,oldsol, t, Δt) end,      # Function for postprocessing successful step
               sample=function(sol,t) end,      # Function to be called for each t\in times[2:end]
               delta=(u,v,t, Δt)->norm(system,u-v,Inf), # Time step error estimator
               transient=true, # choose between transient and stationary (embedding) case
               time=0.0,
               params=zeros(0),
               kwargs...
               )
    
    λstr="t"
    if !transient
        λstr="p"
    end

    allhistory=TransientSolverHistory()
    
    solution=copy(inival)
    oldsolution=copy(inival)
    _initialize_dirichlet!(inival,system;  time, λ=Float64(lambdas[1]),params)
    Δλ=Δλ_val(control,transient)
    
    if !transient
        pre(solution,Float64(lambdas[1]))
        solution=solve!(solution,oldsolution, system; control=control,time=time,tstep=Inf,embedparam=Float64(lambdas[1]),params=params, kwargs...)
        post(solution,oldsolution,lambdas[1], 0)
        if control.log
            push!(allhistory,system.history)
            push!(allhistory.times,lambdas[1])
            Δu=delta(solution, oldsolution,lambdas[1],0)
            push!(allhistory.updates,Δu)
        end
        oldsolution.=solution
    end
    
    tsol=TransientSolution(Float64(lambdas[1]),solution, in_memory=control.in_memory)
    
    if control.verbose
        @printf("  Evolution: start\n")
    end
    
    istep=0
    for i=1:length(lambdas)-1
        
        Δλ=max(Δλ,Δλ_min(control,transient))
        λstart=lambdas[i]
        λend=lambdas[i+1]
        λ=Float64(λstart)
        
        while λ<λend
            solved=false
            λ0=λ
            Δu=0.0
            while !solved
                solved=true
                try
                    λ=λ0+Δλ
                    pre(solution,λ)
                    if transient
                        solution=solve!(solution,oldsolution, system; control=control,time=λ,tstep=Δλ,params=params, kwargs...)
                    else
                        solution=solve!(solution,oldsolution, system; control=control,time=time,tstep=Inf,embedparam=λ,params=params,kwargs...)
                    end
                catch err
                    if (control.handle_exceptions)
                        _print_error(err,stacktrace(catch_backtrace()))
                    else
                        rethrow(err)
                    end
                    solved=false
                end
                if solved
                    Δu=delta(solution, oldsolution,λ, Δλ)
                    if Δu>2.0*control.Δu_opt
                        solved=false
                    end
                end
                if !solved
                    # reduce time step and retry  solution
                    Δλ=Δλ*0.5
                    if Δλ<Δλ_min(control,transient)
                        if !(control.force_first_step && istep==0)
                            throw(EmbeddingError(" Δ$(λstr)_min=$(Δλ_min(control,transient)) reached while Δu=$(Δu) >>  Δu_opt=$(control.Δu_opt) "))
                        else
                            solved=true
                        end
                    end
                    if control.verbose
                        @printf("  Evolution: retry: Δ%s=%.3e\n",λstr,Δλ)
                    end
                end
            end
            istep=istep+1
            if control.verbose
                @printf("  Evolution: istep=%d λ=%.3e Δu=%.3e\n",istep,λ,Δu)
            end
            if control.log
                push!(allhistory,system.history)
                push!(allhistory.updates,Δu)
                push!(allhistory.times,λ)
            end
            if control.store_all
                append!(tsol,λ,solution)
            end
            post(solution,oldsolution,λ, Δλ)
            oldsolution.=solution
            if λ<λend
                Δλ=min(Δλ_max(control,transient),
                       Δλ*Δλ_grow(control,transient),
                       Δλ*control.Δu_opt/(Δu+1.0e-14),
                       λend-λ)
            end
        end
        if !control.store_all
            append!(tsol,lambdas[i+1],solution)
        end
        sample(solution,lambdas[i+1])
    end
    if control.verbose
        @printf("  Evolution: success\n")
    end
    
    system.history=allhistory
    return tsol
end


"""
    evaluate_residual_and_jacobian(system,u;
                                   t=0.0, tstep=Inf,embed=0.0)

Evaluate residual and jacobian at solution value u.
Returns a solution vector containing the residual, and an ExendableSparseMatrix
containing the linearization at u.

"""
function evaluate_residual_and_jacobian(sys,u;t=0.0, tstep=Inf,embed=0.0)
    _complete!(sys, create_newtonvectors=true)

    eval_and_assemble(sys,u,u,sys.residual,t,tstep,embed,zeros(0))
    copy(sys.residual), copy(flush!(sys.matrix))
end


module NoModule end

#####################################################################
"""
    VoronoiFVM.solve(system; kwargs...)
    
Built-in solution method for [`VoronoiFVM.System`](@ref).  
For using ODE solvers from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), see
the [VoronoiFVMDiffEq.jl](https://github.com/j-fu/VoronoiFVMDiffEq.jl) package.
    
Keyword arguments:
- General for all solvers 
   - `inival` (default: 0) : Array created via [`unknowns`](@ref) or  number giving the initial value.
   -  All elements of [`SolverControl`](@ref) can be used as kwargs.
   - `damp` (default: 1): alias for `damp_initial`
   - `damp_grow` (default: 1): alias for `damp_growth`
   - `abstol`: alias for `tol_absolute`
   - `reltol`: alias for `tol_relative`
   - `control` (default: nothing): Pass instance of [`SolverControl`](@ref)
   - `params`: Parameters (Parameter handling is experimental and may change)
    
- __Stationary solver__:
  Invoked if neither `times` nor `embed`, nor `tstep` are given as keyword argument.
  - `time` (default: `0`): Set time value. 
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
  - `delta` (default:  `(u,v,t, Δt)->norm(sys,u-v,Inf)` ):  Value  used to control the time step size `Δu`
  If `control.handle_error` is true, if step solution  throws an error,
  stepsize  is lowered, and  step solution is called again with a smaller time value.
  If `control.Δt<control.Δt_min`, solution is aborted with error.
  Returns a transient solution object `sol` containing the stored solution,  see [`TransientSolution`](@ref).
  
- __Implicit Euler timestep solver__.  Invoked if `tstep` kwarg is given. Solve one time step of the implicit Euler method.
  - `time` (default: `0`): Set time value. 
  - `tstep`: time step
  Returns a [`DenseSolutionArray`](@ref) or [`SparseSolutionArray`](@ref)
"""
function VoronoiFVM.solve(sys::VoronoiFVM.AbstractSystem; inival=0, params=zeros(0), control=VoronoiFVM.NewtonControl(), kwargs...)
    if isa(inival,Number)
        inival=unknowns(sys,inival=inival)
    elseif  !VoronoiFVM.isunknownsof(inival,sys)
        @error "wrong type of inival: $(typeof(inival))"
    end

    # compatibility to names in SolverControl which cannot be deprecated.
    for k ∈ kwargs
        key=k[1]
        if key==:abstol key=:tol_absolute end
        if key==:reltol key=:tol_relative end
        if key==:damp key=:damp_initial end
        if key==:damp_grow key=:damp_growth end
        if hasproperty(control,key)
            setproperty!(control,key,k[2])
        end
    end
    
    if haskey(kwargs,:times)
        solve(inival,sys,kwargs[:times]; control=control, transient=true, params=params,kwargs...)
    elseif haskey(kwargs,:embed)
        solve(inival,sys,kwargs[:embed]; control=control, transient=false,params=params, kwargs...)
    else
        solve(inival,sys; control=control,params=params,kwargs...)
    end
end
