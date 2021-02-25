##################################################################
"""
$(SIGNATURES)

Extract value from dual number. Use to debug physics callbacks.
Re-exported from ForwardDiff.jl
"""
const value=ForwardDiff.value
const partials=ForwardDiff.partials
const npartials=ForwardDiff.npartials



# Add value to matrix if it is nonzero
@inline function _addnz(matrix,i,j,v::Tv,fac) where Tv
    if isnan(v)
        error("trying to assemble NaN")
    end
    if v!=zero(Tv)
        updateindex!(matrix,+,v*fac,i,j)
    end
end



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
end


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

################################################################
function _solve!(
    solution::AbstractMatrix{Tv}, # old time step solution resp. initial value
    oldsol::AbstractMatrix{Tv}, # old time step solution resp. initial value
    this::AbstractSystem{Tv}, # Finite volume system
    control::NewtonControl,
    tstep::Tv,
    log::Bool
) where Tv

    _complete!(this, create_newtonvectors=true)
    nlhistory=zeros(0)
    
    solution.=oldsol
    residual=this.residual
    update=this.update
    _initialize!(solution,this)
    SuiteSparse.UMFPACK.umf_ctrl[3+1]=control.umfpack_pivot_tolerance

    oldnorm=1.0
    converged=false
    if control.verbose
        @printf("    Start Newton iteration\n")
    end
    nlu=0
    nround=0
    lufact=nothing
    damp=control.damp_initial
    tolx=0.0
    rnorm=LinearAlgebra.norm(values(solution),1)

    for ii=1:control.max_iterations
        try
            eval_and_assemble(this,solution,oldsol,residual,tstep,edge_cutoff=control.edge_cutoff)
        catch err
            if (control.handle_exceptions)
                _print_error(err,stacktrace(catch_backtrace()))
                throw(AssemblyError())
            else
                rethrow(err)
            end
        end

        mtx=this.matrix
        
        nliniter=0
        # Sparse LU factorization
        # Here, we seem miss the possibility to re-use the 
        # previous symbolic information
        # We however reuse the factorization control.max_lureuse times.
        if nlu==0
            try
                lufact=LinearAlgebra.lu(mtx)
            catch err
                if (control.handle_exceptions)
                    _print_error(err,stacktrace(catch_backtrace()))
                    throw(FactorizationError())
                else
                    rethrow(err)
                end
            end
            # LU triangular solve gives Newton update
            try
                ldiv!(values(update),lufact,values(residual))
            catch err
                if (control.handle_exceptions)
                    _print_error(err,stacktrace(catch_backtrace()))
                    throw(FactorizationError())
                else
                    rethrow(err)
                end
            end
        else
            # When reusing lu factorization, we may try to iterate
            # Genericly, this is advisable.
            if control.tol_linear <1.0
                (sol,history)= bicgstabl!(values(update),
                                          mtx,
                                          values(residual),
                                          1,
                                          Pl=lufact,
                                          reltol=control.tol_linear,
                                          max_mv_products=100,
                                          log=true)
                nliniter=history.iters
            else
                try
                    ldiv!(values(update),lufact,values(residual))
                catch err
                    if (control.handle_exceptions)
                        _print_error(err,stacktrace(catch_backtrace()))
                        throw(FactorizationError())
                    else
                        rethrow(err)
                    end
                end
            end
        end
        
        nlu=min(nlu+1,control.max_lureuse)
        solval=values(solution)
        solval.-=damp*values(update)
        damp=min(damp*control.damp_growth,1.0)
        norm=LinearAlgebra.norm(values(update),Inf)
        if tolx==0.0
            tolx=norm*control.tol_relative
        end

        dnorm=1.0
        rnorm_new=LinearAlgebra.norm(values(solution),1)
        if rnorm>1.0e-50
            dnorm=abs((rnorm-rnorm_new)/rnorm)
        end
        
        if dnorm<control.tol_round
            nround=nround+1
        else
            nround=0
        end

        push!(nlhistory,norm)
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
    if control.verbose
        @printf("    Newton iteration successful\n")
    end
    return nlhistory
end

################################################################
"""
$(SIGNATURES)

Main assembly method.

Evaluate solution with result in right hand side F and 
assemble matrix into system.matrix.
"""
function eval_and_assemble(system::AbstractSystem{Tv, Ti},
                           U::AbstractMatrix{Tv}, # Actual solution iteration
                           UOld::AbstractMatrix{Tv}, # Old timestep solution
                           F::AbstractMatrix{Tv},# Right hand side
                           tstep::Tv; # time step size. Inf means stationary solution
                           edge_cutoff=0.0
                           ) where {Tv, Ti}
    

    _complete!(system) # needed here as well for test function system which does not use newton
    
    grid=system.grid
    physics=system.physics
    data=physics.data
    node=Node{Tv,Ti}(system)
    bnode=BNode{Tv,Ti}(system)
    edge=Edge{Tv,Ti}(system)
    nspecies=num_species(system)
    matrix=system.matrix
    cellnodefactors=system.cellnodefactors
    celledgefactors=system.celledgefactors
    bfacenodefactors=system.bfacenodefactors
    
    
    # splatting would work but costs allocations
    if isdata(data)
        issource=(physics.source!=nofunc)
        isreaction=(physics.reaction!=nofunc)
        isbreaction=(physics.breaction!=nofunc)
        isbstorage=(physics.bstorage!=nofunc)

        fluxwrap=function(y, u)
            y.=0
            physics.flux(y,u,edge,data)
            nothing
        end

        reactionwrap=function (y, u)
            y.=0
            ## for ii in ..  uu[node.speclist[ii]]=u[ii]
            physics.reaction(y,u,node,data)
            ## for ii in .. y[ii]=y[node.speclist[ii]]
            nothing
        end

        storagewrap= function(y, u)
            y.=0
            physics.storage(y,u,node,data)
            nothing
        end
        
        sourcewrap=function(y)
            y.=0
            physics.source(y,node,data)
            nothing
        end
        
         breactionwrap=function(y, u)
            y.=0
            physics.breaction(y,u,bnode,data)
            nothing
        end
        
        bstoragewrap=function(y, u)
            y.=0
            physics.bstorage(y,u,bnode,data)
            nothing
        end
        
    else
        issource=(physics.source!=nofunc2)
        isreaction=(physics.reaction!=nofunc2)
        isbreaction=(physics.breaction!=nofunc2)
        isbstorage=(physics.bstorage!=nofunc2)

        fluxwrap=function(y, u)
            y.=0
            physics.flux(y,u,edge)
            nothing
        end

        reactionwrap=function(y, u)
            y.=0
            ## for ii in ..  uu[node.speclist[ii]]=u[ii]
            physics.reaction(y,u,node)
            ## for ii in .. y[ii]=y[node.speclist[ii]]
            nothing
        end

        storagewrap= function(y, u)
            y.=0
            physics.storage(y,u,node)
            nothing
        end
        
        sourcewrap=function(y)
            y.=0
            physics.source(y,node)
            nothing
        end
        
         breactionwrap=function(y, u)
            y.=0
            physics.breaction(y,u,bnode)
            nothing
        end
        
        bstoragewrap=function(y, u)
            y.=0
            physics.bstorage(y,u,bnode)
            nothing
        end

    end        
    
    # Reset matrix + rhs
    nzv=nonzeros(matrix)
    nzv.=0.0
    F.=0.0

    # structs holding diff results for storage, reaction,  flux ...
    result_r=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_s=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_flx=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,2*nspecies))
    
    # Array holding function results
    Y=Array{Tv,1}(undef,nspecies)

    # Arrays for gathering solution data
    UK=Array{Tv,1}(undef,nspecies)
    UKOld=Array{Tv,1}(undef,nspecies)
    UKL=Array{Tv,1}(undef,2*nspecies)

    # array holding source term
    src=zeros(Tv,nspecies)

    # arrays holding storage terms for old solution
    oldstor=zeros(Tv,nspecies)
    res_react=zeros(Tv,nspecies)
    jac_react=zeros(Tv,nspecies,nspecies)
    oldbstor=zeros(Tv,nspecies)

    cfg_r=ForwardDiff.JacobianConfig(reactionwrap, Y, UK)
    cfg_s=ForwardDiff.JacobianConfig(storagewrap, Y, UK)
    cfg_br=ForwardDiff.JacobianConfig(breactionwrap, Y, UK)
    cfg_bs=ForwardDiff.JacobianConfig(bstoragewrap, Y, UK)
    cfg_flx=ForwardDiff.JacobianConfig(fluxwrap, Y, UKL)

    
    # Inverse of timestep
    # According to Julia documentation, 1/Inf=0 which
    # comes handy to write compact code here.
    tstepinv=1.0/tstep 

    
    boundary_factors=system.boundary_factors
    boundary_values=system.boundary_values
    
    geom=grid[CellGeometries][1]
    csys=grid[CoordinateSystem]
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    cellregions=grid[CellRegions]
    bgeom=grid[BFaceGeometries][1]
    bfacenodes=grid[BFaceNodes]
    bfaceregions=grid[BFaceRegions]
    nbfaces=num_bfaces(grid)
    ncells=num_cells(grid)
    
    if haskey(grid,CellEdges)
        cellx=grid[CellEdges]
        edgenodes=grid[EdgeNodes]
        has_celledges=true
    else
        cellx=grid[CellNodes]
        edgenodes=local_celledgenodes(geom)
        has_celledges=false
    end

    # Main cell loop
    for icell=1:ncells
        for inode=1:num_nodes(geom)
            _fill!(node,cellnodes,cellregions,inode,icell)
            for ispec=1:nspecies
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
                UK[ispec]=U[ispec,node.index]
                UKOld[ispec]=UOld[ispec,node.index]
            end

            if issource
                sourcewrap(src)
            end
            
            if isreaction
                # Evaluate & differentiate reaction term if present
                ForwardDiff.jacobian!(result_r,reactionwrap,Y,UK,cfg_r)
                res_react=DiffResults.value(result_r)
                jac_react=DiffResults.jacobian(result_r)
            end
            
            # Evaluate & differentiate storage term
            ForwardDiff.jacobian!(result_s,storagewrap,Y,UK,cfg_s)
            res_stor=DiffResults.value(result_s)
            jac_stor=DiffResults.jacobian(result_s)

            # Evaluate storage term for old timestep
            storagewrap(oldstor,UKOld)
            
            
            K=node.index

            for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                ispec=_spec(F,idof,K)
                _add(F,idof,cellnodefactors[inode,icell]*(res_react[ispec]-src[ispec] + (res_stor[ispec]-oldstor[ispec])*tstepinv))
                for jdof=_firstnodedof(F,K):_lastnodedof(F,K)
                    jspec=_spec(F,jdof,K)
                    _addnz(matrix,idof,jdof,jac_react[ispec,jspec]+ jac_stor[ispec,jspec]*tstepinv,cellnodefactors[inode,icell])
                end
            end
        end
        for iedge=1:num_edges(geom)
            if abs(celledgefactors[iedge,icell])<edge_cutoff
                continue
            end

            _fill!(edge,cellx,edgenodes,cellregions,iedge,icell, has_celledges)

            #Set up argument for fluxwrap
            for ispec=1:nspecies
                UKL[ispec]=U[ispec,edge.node[1]]
                UKL[nspecies+ispec]=U[ispec,edge.node[2]]
            end

            ForwardDiff.jacobian!(result_flx,fluxwrap,Y,UKL,cfg_flx)
            res=DiffResults.value(result_flx)
            jac=DiffResults.jacobian(result_flx)

            K=edge.node[1]
            L=edge.node[2]
            fac=celledgefactors[iedge,icell]
            for idofK=_firstnodedof(F,K):_lastnodedof(F,K)
                ispec=_spec(F,idofK,K)
                idofL=dof(F,ispec,L)
                if idofL==0
                    continue
                end
                _add(F,idofK,fac*res[ispec])
                _add(F,idofL,-fac*res[ispec])
                
                for jdofK=_firstnodedof(F,K):_lastnodedof(F,K)
                    jspec=_spec(F,jdofK,K)
                    jdofL=dof(F,jspec,L)
                    if jdofL==0
                        continue
                    end
                    
                    _addnz(matrix,idofK,jdofK,+jac[ispec,jspec         ],fac)
                    _addnz(matrix,idofL,jdofK,-jac[ispec,jspec         ],fac)
                    _addnz(matrix,idofK,jdofL,+jac[ispec,jspec+nspecies],fac)
                    _addnz(matrix,idofL,jdofL,-jac[ispec,jspec+nspecies],fac)
                    
                end
            end
        end
    end

    # Assembly loop for boundary conditions
    for ibface=1:nbfaces
        # Obtain boundary region number
        ibreg=bfaceregions[ibface]

        # Loop over nodes of boundary face
        for ibnode=1:num_nodes(bgeom)

            # Fill bnode data shuttle with data from grid
            _fill!(bnode,bfacenodes,bfaceregions,ibnode,ibface)

            # Copy unknown values from solution into dense array
            for ispec=1:nspecies
                UK[ispec]=U[ispec,bnode.index]
            end

            # Measure of boundary face part assembled to node
            bnode_factor=bfacenodefactors[ibnode,ibface]

            # Global index of node
            K=bnode.index
            
            # Assemble "standard" boundary conditions: Robin or
            # Dirichlet
            # valid only for interior species, currently not checked
            for ispec=1:nspecies
                idof=dof(F,ispec,K)
                # If species is present, assemble the boundary condition
                if idof>0
                    # Get user specified data
                    boundary_factor=boundary_factors[ispec,ibreg]
                    boundary_value=boundary_values[ispec,ibreg]
                    
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

            # Boundary reaction term has been given.
            # Valid for both boundary and interior species
            if isbreaction

                # Evaluate function + derivative
                ForwardDiff.jacobian!(result_r,breactionwrap,Y,UK,cfg_br)
                res_breact=DiffResults.value(result_r)
                jac_breact=DiffResults.jacobian(result_r)
                
                # Assemble RHS + matrix
                for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                    ispec=_spec(F,idof,K)
                    _add(F,idof,bnode_factor*res_breact[ispec])
                    for jdof=_firstnodedof(F,K):_lastnodedof(F,K)
                        jspec=_spec(F,jdof,K)
                        _addnz(matrix,idof,jdof,jac_breact[ispec,jspec],bnode_factor)
                    end
                end

            end
            
            # Boundary reaction term has been given.
            # Valid only for boundary species, but this is currently not checked.
            if isbstorage
                
                # Fetch data for old timestep
                for ispec=1:nspecies
                    UKOld[ispec]=UOld[ispec,K]
                end
                # Evaluate storage term for old timestep
                bstoragewrap(oldbstor,UKOld)

                # Evaluate & differentiate storage term for new time step value
                ForwardDiff.jacobian!(result_s,bstoragewrap,Y,UK,cfg_bs)
                res_bstor=DiffResults.value(result_s)
                jac_bstor=DiffResults.jacobian(result_s)
                
                for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                    ispec=_spec(F,idof,K)
                    # Assemble finite difference in time for right hand side
                    _add(F,idof,bnode_factor*(res_bstor[ispec]-oldbstor[ispec])*tstepinv)
                    # Assemble matrix.
                    for jdof=_firstnodedof(F,K):_lastnodedof(F,K)
                        jspec=_spec(F,jdof,K)
                        _addnz(matrix,idof,jdof,jac_bstor[ispec,jspec],bnode_factor*tstepinv)
                    end
                end
            end
            
        end
    end
    _eval_and_assemble_generic_operator(system,U,F)
    _eval_and_assemble_inactive_species(system,U,UOld,F)
end


function _eval_and_assemble_generic_operator(this::AbstractSystem,U,F)

    if !has_generic_operator(this)
        return
    end
    generic_operator(f,u)=this.physics.generic_operator(f,u,this)
    vecF=values(F)
    vecU=values(U)
    y=similar(vecF)
    generic_operator(y,vecU)
    vecF.+=y
    forwarddiff_color_jacobian!(this.generic_matrix, generic_operator, vecU, colorvec = this.generic_matrix_colors)
    rowval=this.generic_matrix.rowval
    colptr=this.generic_matrix.colptr
    nzval=this.generic_matrix.nzval
    for i=1:length(colptr)-1
        for j=colptr[i]:colptr[i+1]-1
            updateindex!(this.matrix,+,nzval[j],i,rowval[j])
        end
    end
end

################################################################
"""
$(TYPEDSIGNATURES)

Solution method for instance of abstract system.

Perform solution of stationary system (if `tstep==Inf`) or one tine step
of implicit Euler time step system using Newton's method with damping. 
Initial damping is chosen according  `control.damp_initial`.
"""
function solve!(
    solution::AbstractMatrix{Tv}, # Solution
    inival::AbstractMatrix{Tv},   # Initial value 
    this::AbstractSystem{Tv};     # Finite volume system
    control=NewtonControl(),      # Newton solver control information
    tstep::Tv=Inf,                # Time step size. Inf means  stationary solution
    log::Bool=false
) where Tv
    if control.verbose
        @time begin
            history=_solve!(solution,inival,this,control,tstep,log)
        end
    else 
        history=_solve!(solution,inival,this,control,tstep,log)
    end
    return log ? (solution,history) : solution
end


################################################################
"""
$(TYPEDSIGNATURES)

Solution method for instance of abstract system.

Perform solution of stationary system (if `tstep==Inf`) or one time step
of implicit Euler time step system using Newton's method with damping. 
Initial damping is chosen according  `control.damp_initial`.
"""
function solve(
    inival::AbstractMatrix{Tv},   # Initial value 
    this::AbstractSystem{Tv};     # Finite volume system
    control=NewtonControl(),      # Newton solver control information
    tstep::Tv=Inf,                # Time step size. Inf means  stationary solution
    log::Bool=false
) where Tv
    solve!(unknowns(this),inival,this,control=control, tstep=tstep,log=log)
end

################################################################
"""
$(SIGNATURES)

Solution method for instance of abstract system.

Perform solution via parameter embedding, calling
solve! for each value of the parameter p from interval
(0,1). The user is responsible for the interpretation of
the parameter. The optional `pre()` callback can be used
to communicate its value to the system.
The optional `post()` callback method can be used to perform
various postprocessing steps.

If `control.handle_error` is true, `solve!`  throws an error, and
 stepsize `control.Δp` is lowered,
and `solve!` is called again with a smaller  parameter
value. If `control.Δp<control.Δp_min`, `embed!` is aborted
with error.

"""
function embed!(
    solution::AbstractMatrix{Tv}, # Solution
    xinival::AbstractMatrix{Tv},   # Initial value 
    this::AbstractSystem{Tv};     # Finite volume system
    control=NewtonControl(),      # Newton solver control information
    pre=function(sol,p) end,       # Function for preparing step
    post=function(sol,p) end      # Function for postprocessing successful step
) where Tv
    inival=copy(xinival)
    p=0.0
    Δp=control.Δp
    if control.verbose
        @printf("  Embedding: start\n")
    end
    istep=0
    while p<1.0
        solved=false
        p0=p
        while !solved
            solved=true
            try
                p=min(p0+Δp,1.0)
                pre(solution,p)
                solve!(solution,inival, this ,control=control)
            catch err
                if (control.handle_exceptions)
                    _print_error(err,stacktrace(catch_backtrace()))
                else
                    rethrow(err)
                end
                # here we reduce the embedding step and retry the solution
                Δp=Δp*0.5
                if Δp<control.Δp_min
                    throw(EmbeddingError())
                end
                solved=false
                if control.verbose
                    @printf("  Embedding: retry: Δp=%.3e\n",Δp)
                end
            end
        end
        istep=istep+1
        if control.verbose
            @printf("  Embedding: istep=%d p=%.3e\n",istep,p)
        end
        Δp=min(control.Δp_max,Δp*1.2)
        post(solution,p)
        inival.=solution
    end
    if control.verbose
        @printf("  Embedding: success\n")
    end
    return solution
end

################################################################
"""
$(TYPEDSIGNATURES)

Time dependent solver for abstract system

Use implicit Euler method with + damped   Newton's method  to 
solve time dependent problem.
"""
function evolve!(
    solution::AbstractMatrix{Tv}, # Solution
    xinival::AbstractMatrix{Tv},   # Initial value 
    this::AbstractSystem{Tv},    # Finite volume system
    times::AbstractVector;
    control=NewtonControl(),      # Newton solver control information
    pre=function(sol,t) end,       # Function for preparing step
    post=function(sol,oldsol, t, Δt) end,      # Function for postprocessing successful step
    delta=(u,v,t, Δt)->norm(u-v,Inf) # Time step error estimator
) where Tv
    inival=copy(xinival)
    Δt=control.Δt
    if control.verbose
        @printf("  Evolution: start\n")
    end
    for i=1:length(times)-1
        Δt=max(Δt,control.Δt_min)
        tstart=times[i]
        tend=times[i+1]
        t=tstart
        istep=0

        while t<tend
            solved=false
            t0=t
            Δu=0.0
            while !solved
                solved=true
                try
                    t=t0+Δt
                    pre(solution,t)
                    solve!(solution,inival, this ,control=control,tstep=Δt)
                catch err
                    if (control.handle_exceptions)
                        _print_error(err,stacktrace(catch_backtrace()))
                    else
                        rethrow(err)
                    end
                    solved=false
                end
                if solved
                    Δu=delta(solution, inival,t, Δt)
                    if Δu>2.0*control.Δu_opt
                        solved=false
                    end
                end
                if !solved
                    # reduce time step and retry  solution
                    Δt=Δt*0.5
                    if Δt<control.Δt_min
                        @printf(" Δt_min=%.2g reached while Δu=%.2g >>  Δu_opt=%.2g\n",control.Δt_min, Δu,control.Δu_opt)
                        throw(EmbeddingError())
                    end
                    if control.verbose
                        @printf("  Evolution: retry: Δt=%.3e\n",Δt)
                    end
                end
            end
            istep=istep+1
            if control.verbose
                @printf("  Evolution: istep=%d t=%.3e Δu=%.3e\n",istep,t,Δu)
            end
            post(solution,inival,t, Δt)
            inival.=solution
            if t<tend
                Δt=min(control.Δt_max,
                       Δt*control.Δt_grow,
                       Δt*control.Δu_opt/(Δu+1.0e-14),
                       tend-t)
            end
        end
    end
    if control.verbose
        @printf("  Evolution: success\n")
    end
    return solution
end



################################################################
"""
$(SIGNATURES)

Integrate function `F` of  solution vector over domain or boundary 
The result contains the integral for each species separately.
"""
function integrate(this::AbstractSystem{Tv,Ti,Tm},F::Function,U::AbstractMatrix{Tu}; boundary=false) where {Tu,Tv,Ti,Tm}
    grid=this.grid
    data=this.physics.data
    nspecies=num_species(this)
    res=zeros(Tu,nspecies)
    node=Node{Tv,Ti}(this)
    nodeparams=(node,)
    if isdata(data)
        nodeparams=(node,data,)
    end

    

    csys=grid[CoordinateSystem]
    coord=grid[Coordinates]


    if boundary
        
        geom=grid[BFaceGeometries][1]
        bfacenodes=grid[BFaceNodes]
        bfaceregions=grid[BFaceRegions]
        nbfaceregions=maximum(bfaceregions)
        integral=zeros(Tu,nspecies,nbfaceregions)
        
        
        for ibface=1:num_bfaces(grid)
            for inode=1:num_nodes(geom)
                _fill!(node,bfacenodes,bfaceregions,inode,ibface)
                F(res,U[:,node.index],nodeparams...)
                for ispec=1:nspecies
                    if this.node_dof[ispec,node.index]==ispec
                        integral[ispec,node.region]+this.bfacenodefactors[inode,ibface]*res[ispec]
                    end
                end
            end
        end
    else
        geom=grid[CellGeometries][1]
        cellnodes=grid[CellNodes]
        cellregions=grid[CellRegions]
        ncellregions=maximum(cellregions)
        integral=zeros(Tu,nspecies,ncellregions)
        
        
        for icell=1:num_cells(grid)
            for inode=1:num_nodes(geom)
                _fill!(node,cellnodes,cellregions,inode,icell)
                F(res,U[:,node.index],nodeparams...)
                for ispec=1:nspecies
                    if this.node_dof[ispec,node.index]==ispec
                        integral[ispec,node.region]+=this.cellnodefactors[inode,icell]*res[ispec]
                    end
                end
            end
        end
    end
    

    return integral
end

