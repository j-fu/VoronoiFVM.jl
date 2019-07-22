# Add value to matrix if it is nonzero
@inline function _addnz(matrix,i,j,v::Tv,fac) where Tv
    if v!=zero(Tv)
        matrix[i,j]+=v*fac
    end
end

# Main assembly loop
function _eval_and_assemble(this::AbstractSystem{Tv},
                            U::AbstractMatrix{Tv}, # Actual solution iteration
                            UOld::AbstractMatrix{Tv}, # Old timestep solution
                            F::AbstractMatrix{Tv},# Right hand side
                            tstep::Tv, # time step size. Inf means stationary solution
                            xstorage::FSTOR, # ensure type stability
                            xbstorage::FBSTOR, # ensure type stability
                            xsource::FSRC,  # ensure type stability
                           ) where {Tv,FSRC,FSTOR, FBSTOR}

    _complete!(this) # needed here as well for test function system which does not use newton
    
    grid=this.grid

    physics::Physics=this.physics
    data=physics.data
    node::Node{Tv}=Node{Tv}(this)
    bnode::BNode{Tv}=BNode{Tv}(this)
    edge::Edge{Tv}=Edge{Tv}(this)
    edge_cutoff=1.0e-12
    nspecies::Int32=num_species(this)
    
    matrix=this.matrix
    issource=(physics.source!=nofunc)
    isreaction=(physics.reaction!=nofunc)
    isbreaction=(physics.breaction!=nofunc)
    isbstorage=(physics.bstorage!=nofunc)

    
    function fluxwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.flux(y,u,edge,data)
    end


    function reactionwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        ## for ii in ..  uu[node.speclist[ii]]=u[ii]
        physics.reaction(y,u,node,data)
        ## for ii in .. y[ii]=y[node.speclist[ii]]
    end

    function storagewrap(y::AbstractVector, u::AbstractVector)
        y.=0
        xstorage(y,u,node,data)
    end

    function sourcewrap(y)
        y.=0
        xsource(y,node,data)
    end

    @inline function breactionwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.breaction(y,u,bnode,data)
    end

    @inline function bstoragewrap(y::AbstractVector, u::AbstractVector)
        y.=0
        xbstorage(y,u,bnode,data)
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

    
    # Arrays holding for factor data
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))
    bnode_factors=zeros(Tv,num_nodes_per_bface(grid))

    # Main cell loop
    for icell=1:num_cells(grid)
        # set up form factors
        cellfactors!(grid,icell,node_factors,edge_factors)
        # set up data for callbacks

        node.region=reg_cell(grid,icell)
        edge.region=reg_cell(grid,icell)

        for inode=1:num_nodes_per_cell(grid)
            _fill!(node,grid, inode,icell)
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
                # Evaluate reaction term if present
                ForwardDiff.jacobian!(result_r,reactionwrap,Y,UK,cfg_r)
                res_react=DiffResults.value(result_r)
                jac_react=DiffResults.jacobian(result_r)
            end
            # Evaluate & differentiate storage term
            ForwardDiff.jacobian!(result_s,storagewrap,Y,UK,cfg_s)
            res_stor=DiffResults.value(result_s)
            jac_stor=DiffResults.jacobian(result_s)

            # Evaluate storage term for old timestep
            ### ALLOC ???
            storagewrap(oldstor,UKOld) 
            
            
            K=node.index

            for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                ispec=_spec(F,idof,K)
                _add(F,idof,node_factors[inode]*(res_react[ispec]-src[ispec] + (res_stor[ispec]-oldstor[ispec])*tstepinv))
                for jdof=_firstnodedof(F,K):_lastnodedof(F,K)
                    jspec=_spec(F,jdof,K)
                    _addnz(matrix,idof,jdof,jac_react[ispec,jspec]+ jac_stor[ispec,jspec]*tstepinv,node_factors[inode])
                end
            end
        end
        
        for iedge=1:num_edges_per_cell(grid)
            if edge_factors[iedge]<edge_cutoff
                continue
            end

            _fill!(edge,grid,iedge,icell)

            #Set up argument for fluxwrap
            for ispec=1:nspecies
                UKL[ispec]=U[ispec,edge.nodeK]
                UKL[ispec+nspecies]=U[ispec,edge.nodeL]
            end

            ForwardDiff.jacobian!(result_flx,fluxwrap,Y,UKL,cfg_flx)
            res=DiffResults.value(result_flx)
            jac=DiffResults.jacobian(result_flx)

            K=edge.nodeK
            L=edge.nodeL
            fac=edge_factors[iedge]
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
                    
                    _addnz(matrix,idofK,jdofK,+jac[ispec,jspec            ],fac)
                    _addnz(matrix,idofL,jdofK,-jac[ispec,jspec            ],fac)
                    _addnz(matrix,idofK,jdofL,+jac[ispec,jspec+nspecies],fac)
                    _addnz(matrix,idofL,jdofL,-jac[ispec,jspec+nspecies],fac)
                    
                end
            end
        end
    end

    for ibface=1:num_bfaces(grid)
        bfacefactors!(grid,ibface,bnode_factors)
        ibreg=grid.bfaceregions[ibface]

        for ibnode=1:num_nodes_per_bface(grid)
            _fill!(bnode,grid,ibnode,ibface)
            for ispec=1:nspecies
                UK[ispec]=U[ispec,bnode.index]
                UKOld[ispec]=UOld[ispec,bnode.index]
            end
            for ispec=1:nspecies # should involve only rspecies
                idof=dof(F,ispec,bnode.index)
                if idof>0
                    fac=this.boundary_factors[ispec,ibreg]
                    val=this.boundary_values[ispec,ibreg]
                    if fac==Dirichlet
                        F[ispec,bnode.index]+=fac*(U[ispec,bnode.index]-val)
                        _addnz(matrix,idof,idof,fac,1)
                    else
                        F[ispec,bnode.index]+=bnode_factors[ibnode]*(fac*U[ispec,bnode.index]-val)
                        _addnz(matrix,idof,idof,fac,bnode_factors[ibnode])
                    end
                end
            end

            if isbreaction# involves bspecies and species
                ForwardDiff.jacobian!(result_r,breactionwrap,Y,UK,cfg_br)
                res_breact=DiffResults.value(result_r)
                jac_breact=DiffResults.jacobian(result_r)
                K=bnode.index
                fac=bnode_factors[ibnode]
                for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                    ispec=_spec(F,idof,K)
                    _add(F,idof,fac*res_breact[ispec])
                    for jdof=_firstnodedof(F,K):_lastnodedof(F,K)
                        jspec=_spec(F,jdof,K)
                        _addnz(matrix,idof,jdof,jac_breact[ispec,jspec],fac)
                    end
                end

            end
            
            if isbstorage # should involve only bspecies
                # Evaluate & differentiate storage term
                ForwardDiff.jacobian!(result_s,bstoragewrap,Y,UK,cfg_bs)
                res_bstor=DiffResults.value(result_s)
                jac_bstor=DiffResults.jacobian(result_s)
                
                # Evaluate storage term for old timestep
                bstoragewrap(oldbstor,UKOld)
                
                for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                    ispec=_spec(F,idof,K)
                    _add(F,idof,fac*(res_bstor[ispec]-oldbstor[ispec])*tstepinv)
                    for jdof=_firstnodedof(F,K):_lastnodedof(F,K)
                        jspec=_spec(F,jdof,K)
                        _addnz(matrix,idof,jdof,jac_bstor[ispec,jspec],fac*tstepinv)
                    end
                end
            end
            
        end
    end
    _eval_and_assemble_inactive_species(this,U,UOld,F)
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
    tstep::Tv
) where Tv

    _complete!(this, create_newtonvectors=true)

    solution.=oldsol
    residual=this.residual
    update=this.update
    _initialize!(solution,this)

    # Newton iteration
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
            _eval_and_assemble(this,solution,oldsol,residual,tstep,
                               this.physics.storage,
                               this.physics.bstorage,
                               this.physics.source)
        catch err
            if (control.handle_exceptions)
                _print_error(err,stacktrace(catch_backtrace()))
                throw(AssemblyError())
            else
                rethrow(err)
            end
        end

        
        # Sparse LU factorization
        # Here, we seem miss the possibility to re-use the 
        # previous symbolic information
        # We however reuse the factorization control.max_lureuse times.
        nliniter=0
        if nlu==0
            lufact=LinearAlgebra.lu(this.matrix)
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
            # Generally, this is advisable.
            if control.tol_linear <1.0
                (sol,history)= bicgstabl!(values(update),
                                          this.matrix,
                                          values(residual),
                                          1,
                                          Pl=lufact,
                                          tol=control.tol_linear,
                                          max_mv_products=20,
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
end

################################################################
"""
$(TYPEDSIGNATURES)

Wrapper for main assembly method.

Evaluate solution with result in right hand side F and 
assemble matrix into system.matrix.
"""
function eval_and_assemble(system::AbstractSystem{Tv},
                           U::AbstractMatrix{Tv}, # Actual solution iteration
                           UOld::AbstractMatrix{Tv}, # Old timestep solution
                           F::AbstractMatrix{Tv},# Right hand side
                           tstep::Tv, # time step size. Inf means stationary solution
                           ) where {Tv}

    _eval_and_assemble(system,
                       solution,
                       oldsol,
                       residual,
                       tstep,
                       system.physics.storage,
                       system.physics.bstorage,
                       system.physics.source)
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
    tstep::Tv=Inf                 # Time step size. Inf means  stationary solution
) where Tv
    if control.verbose
        @time begin
            _solve!(solution,inival,this,control,tstep)
        end
    else 
        _solve!(solution,inival,this,control,tstep)
    end
    return solution
end



################################################################
"""
$(TYPEDSIGNATURES)

Solution method for instance of abstract system.

Perform solution via parameter embedding, calling
solve! for each value of the parameter p from interval
(0,1). The user is responsible for the interpretation of
the parameter. The optional `pre()` callback can be used
to communicate its value to the system.
The optional`post()` callback method can be used to perform
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

Integrate function `F` of  solution vector over domain. 
The result contains the integral for each species separately.
"""
function integrate(this::AbstractSystem{Tv},F::Function,U::AbstractMatrix{Tv})::Array{Tv,1} where Tv
    grid=this.grid
    nspecies=num_species(this)
    integral=zeros(Tv,nspecies)
    res=zeros(Tv,nspecies)
    node=Node{Tv}(this)
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))

    for icell=1:num_cells(grid)
        cellfactors!(grid,icell,node_factors,edge_factors)
        for inode=1:num_nodes_per_cell(grid)
            _fill!(node,grid,inode,icell)
            F(res,U[:,node.index],node,this.physics.data)
            for ispec=1:nspecies
                if this.node_dof[ispec,node.index]==ispec
                    integral[ispec]+=node_factors[inode]*res[ispec]
                end
            end
        end
    end
    return integral
end

