##################################################################
"""
$(SIGNATURES)

Extract value from dual number. Use to debug physics callbacks.
Re-exported from ForwardDiff.jl
"""
const value=ForwardDiff.value

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
for implicit Euler and stationary solve.
"""
function _solve!(
    solution::AbstractMatrix{Tv}, # old time step solution resp. initial value
    oldsol::AbstractMatrix{Tv}, # old time step solution resp. initial value
    system::AbstractSystem{Tv}, # Finite volume system
    control::NewtonControl,
    tstep::Tv,
    log::Bool
) where Tv

    _complete!(system, create_newtonvectors=true)
    nlhistory=zeros(0)
    
    solution.=oldsol
    residual=system.residual
    update=system.update
    _initialize!(solution,system)
    SuiteSparse.UMFPACK.umf_ctrl[3+1]=control.umfpack_pivot_tolerance

    oldnorm=1.0
    converged=false
    if control.verbose
        @printf("    Start Newton iteration\n")
    end
    nlu=0
    nround=0
    damp=control.damp_initial
    tolx=0.0
    rnorm=LinearAlgebra.norm(values(solution),1)

    for ii=1:control.max_iterations
        try
            eval_and_assemble(system,solution,oldsol,residual,tstep,edge_cutoff=control.edge_cutoff)
        catch err
            if (control.handle_exceptions)
                _print_error(err,stacktrace(catch_backtrace()))
                throw(AssemblyError())
            else
                rethrow(err)
            end
        end

        mtx=system.matrix
        
        nliniter=0
        # Sparse LU factorization
        if nlu==0 # (re)factorize, if possible reusing old factorization data
            try
                system.factorization=factorize!(system.factorization,mtx;kind=control.factorization)
            catch err
                if (control.handle_exceptions)
                    _print_error(err,stacktrace(catch_backtrace()))
                    throw(FactorizationError())
                else
                    rethrow(err)
                end
            end
        end
        if issolver(system.factorization) && nlu==0
            # Direct solution via LU solve
            try
                ldiv!(values(update),system.factorization,values(residual))
            catch err
                if (control.handle_exceptions)
                    _print_error(err,stacktrace(catch_backtrace()))
                    throw(FactorizationError())
                else
                    rethrow(err)
                end
            end
        elseif !issolver(system.factorization) || nlu>0
            # Iterative solution
            update.=zero(Tv)
            (sol,history)= bicgstabl!(values(update),
                                      mtx,
                                      values(residual),
                                      1,
                                      Pl=system.factorization,
                                      reltol=control.tol_linear,
                                      max_mv_products=100,
                                      log=true)
            nliniter=history.iters
        else
            error("This should not happen")
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
assemble Jacobi matrix into system.matrix.
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
    node=Node{Tv,Ti}(system)
    bnode=BNode{Tv,Ti}(system)
    edge=Edge{Tv,Ti}(system)
    @create_physics_wrappers(physics,node,bnode,edge)


    nspecies::Int=num_species(system)
    matrix=system.matrix
    cellnodefactors::Array{Tv,2}=system.cellnodefactors
    celledgefactors::Array{Tv,2}=system.celledgefactors
    bfacenodefactors::Array{Tv,2}=system.bfacenodefactors
    
    
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

    
    boundary_factors::Array{Tv,2}=system.boundary_factors
    boundary_values::Array{Tv,2}=system.boundary_values
    
    geom=grid[CellGeometries][1]
    csys=grid[CoordinateSystem]
    coord=grid[Coordinates]
    cellnodes::Array{Ti,2}=grid[CellNodes]
    cellregions::Vector{Ti}=grid[CellRegions]
    bgeom=grid[BFaceGeometries][1]
    bfacenodes::Array{Ti,2}=grid[BFaceNodes]
    bfaceregions::Vector{Ti}=grid[BFaceRegions]
    nbfaces=num_bfaces(grid)
    ncells=num_cells(grid)


    cellx::Array{Ti,2}=grid[CellNodes]
    edgenodes::Array{Ti,2}=local_celledgenodes(geom)
    has_celledges=false
    if haskey(grid,CellEdges)
        cellx=grid[CellEdges]
        edgenodes=grid[EdgeNodes]
        has_celledges=true
    end

    nn::Int=num_nodes(geom)
    ne::Int=num_edges(geom)


    ncalloc=@allocated  for icell=1:ncells
        for inode=1:nn
            node.region=cellregions[icell]
            node.index=cellnodes[inode,icell]
            node.icell=icell
            @views UK.=U[:,node.index]
            @views UKOld.=UOld[:,node.index]
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


           if issource
                sourcewrap(src)
            end
            
            
           if isreaction
                # Evaluate & differentiate reaction term if present
                Y.=zero(Tv)
                ForwardDiff.vector_mode_jacobian!(result_r,reactionwrap,Y,UK,cfg_r)
                res_react=DiffResults.value(result_r)
                jac_react=DiffResults.jacobian(result_r)
            end
            
            # Evaluate & differentiate storage term
            Y.=zero(Tv)
            ForwardDiff.vector_mode_jacobian!(result_s,storagewrap,Y,UK,cfg_s)
            res_stor=DiffResults.value(result_s)
            jac_stor=DiffResults.jacobian(result_s)

            # Evaluate storage term for old timestep
            oldstor.=zero(Tv)
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

        for iedge=1:ne
            if abs(celledgefactors[iedge,icell])<edge_cutoff
                continue
            end
            
            if has_celledges #  cellx==celledges, edgenodes==global_edgenodes
                # If we work with projections of fluxes onto edges,
                # we need to ensure that the edges are accessed with the
                # same orientation without regard of the orientation induced
                # by local cell numbering
                edge.index=cellx[iedge,icell]
                edge.node[1]=edgenodes[1,edge.index]
                edge.node[2]=edgenodes[2,edge.index]
            else # cx==cellnodes, edgenodes== local_edgenodes
                edge.index=0
                edge.node[1]=cellx[edgenodes[1,iedge],icell]
                edge.node[2]=cellx[edgenodes[2,iedge],icell]
            end
            edge.region=cellregions[icell]
            edge.icell=icell
            

            #Set up argument for fluxwrap
            @views UKL[1:nspecies].=U[:,edge.node[1]]
            @views UKL[nspecies+1:2*nspecies].=U[:,edge.node[2]]

            Y.=zero(Tv)
            ForwardDiff.vector_mode_jacobian!(result_flx,fluxwrap,Y,UKL,cfg_flx)
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

    nbn::Int=num_nodes(bgeom)
    # Assembly loop for boundary conditions
    nballoc= @allocated for ibface=1:nbfaces
        ibreg::Int=bfaceregions[ibface]

        # Loop over nodes of boundary face
        for ibnode=1:nbn
            # Fill bnode data shuttle with data from grid
            bnode.ibface=ibface
            bnode.ibnode=ibnode
            bnode.region=bfaceregions[ibface]
            bnode.index=bfacenodes[ibnode,ibface]

            # Copy unknown values from solution into dense array
            @views UK.=U[:,bnode.index]

            # Measure of boundary face part assembled to node
            bnode_factor::Tv=bfacenodefactors[ibnode,ibface]

            # Global index of node
            K::Int=bnode.index
            
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
                Y.=zero(Tv)
                ForwardDiff.vector_mode_jacobian!(result_r,breactionwrap,Y,UK,cfg_br)
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
            # Valid only for boundary species, but system is currently not checked.
            if isbstorage
                
                # Fetch data for old timestep
                @views UKOld.=UOld[:,bnode.index]

                # Evaluate storage term for old timestep
                bstoragewrap(oldbstor,UKOld)

                # Evaluate & differentiate storage term for new time step value
                Y.=zero(Tv)
                ForwardDiff.vector_mode_jacobian!(result_s,bstoragewrap,Y,UK,cfg_bs)
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

   if isnothing(matrix.lnkmatrix)
       # There should be no allocations if the matrix pattern has not
       # changed - this is the case when no new entries have been
       # collected into lnkmatrix.
       if system.allocs>=0 # we had a couple of runs before to bridge the compilation phase
           system.allocs=ncalloc+nballoc
           if system.allocs>0
               error("""Allocations in assembly loop: cells: $(ncalloc), bfaces: $(nballoc)
                        See the documentation of `check_allocs!` for more information""")
           end
       elseif system.allocs > -100 # probably still in compiling phase
           system.allocs=system.allocs+1
       end
       # otherwise, checking has been switched off.
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
            updateindex!(system.matrix,+,nzval[j],i,rowval[j])
        end
    end
end

################################################################
"""
````
solve!(solution, inival, system; 
    control=NewtonControl(), 
    tstep=Inf, 
    log=false)
````

Perform solution of stationary problem(if `tstep==Inf`) or of  one step
of the implicit Euler method using Newton's method with `inival` as initial
value. The method writes into `solution`. 

It returns `solution` or, if `log==true`, a tuple of solution and a vector
containing the residual history of Newton's method.
"""
function solve!(
    solution::AbstractMatrix{Tv}, # Solution
    inival::AbstractMatrix{Tv},   # Initial value 
    system::AbstractSystem{Tv};     # Finite volume system
    control=NewtonControl(),      # Newton solver control information
    tstep::Tv=Inf,                # Time step size. Inf means  stationary solution
    log::Bool=false
) where Tv
    if control.verbose
        @time begin
            history=_solve!(solution,inival,system,control,tstep,log)
        end
    else 
        history=_solve!(solution,inival,system,control,tstep,log)
    end
    return log ? (solution,history) : solution
end


################################################################
"""
````
solve(inival, system; 
      control=NewtonControl(), 
      tstep=Inf, 
      log=false)
````

Perform solution of stationary problem(if `tstep==Inf`) or one step
of the implicit Euler method using Newton's method with `inival` as initial
value.
It returns a solution array or, if `log==true`, a tuple of solution and a vector
containing the residual history of Newton's method.
"""
function solve(
    inival::AbstractMatrix{Tv},   # Initial value 
    system::AbstractSystem{Tv};     # Finite volume system
    control=NewtonControl(),      # Newton solver control information
    tstep::Tv=Inf,                # Time step size. Inf means  stationary solution
    log::Bool=false
) where Tv
    solve!(unknowns(system),inival,system,control=control, tstep=tstep,log=log)
end

################################################################
"""
````
function embed!(solution, inival, system; 
                control=NewtonControl(),
                pre=function(sol,p) end,
                post=function(sol,p) end
)
````
Solve stationary problem via parameter embedding, calling
`solve!` for increasing values of the parameter p from interval
(0,1). The user is responsible for the interpretation of
the parameter. 

The optional `pre()` callback can be used
to communicate its value to the system.

The optional `post()` callback method can be used to perform
various postprocessing steps.

If `control.handle_error` is true, if `solve!`  throws an error,
stepsize `control.Δp` is lowered, and `solve!` is called again with a smaller  parameter
value. If `control.Δp<control.Δp_min`, `embed!` is aborted
with error.

"""
function embed!(
    solution::AbstractMatrix{Tv}, # Solution
    xinival::AbstractMatrix{Tv},   # Initial value 
    system::AbstractSystem{Tv};     # Finite volume system
    control=NewtonControl(),      # Newton solver control information
    pre=function(sol,p) end,       # Function for preparing step
    post=function(sol,p) end      # Function for postprocessing successful step
) where Tv
    inival=copy(inival)
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
                solve!(solution,inival, system ,control=control)
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
````
function evolve!(solution, inival, system, times;
                 control=NewtonControl(), 
                 pre=function(sol,t) end,   
                 post=function(sol,oldsol, t, Δt) end,
                 sample=function(sol,t) end,
                 delta=(u,v,t, Δt)->norm(sys,u-v,Inf)
)
````

Use implicit Euler method  + damped   Newton's method  to 
solve time dependent problem. Time step control is performed
according to the data in `control`.  All times in `times`
are reached exactly.

Callbacks:
- `pre` is invoked before each time step
- `post`  is invoked after each time step
- `sample` is called for all times in `times[2:end]`.

`delta` is  used to control the time step.
"""
function evolve!(
    solution::AbstractMatrix{Tv}, # Solution
    inival::AbstractMatrix{Tv},   # Initial value 
    system::AbstractSystem{Tv},    # Finite volume system
    times::AbstractVector;
    control=NewtonControl(),      # Newton solver control information
    pre=function(sol,t) end,       # Function for preparing step
    post=function(sol,oldsol, t, Δt) end,      # Function for postprocessing successful step
    sample=function(sol,t) end,      # Function to be called for each t\in times[2:end]
    delta=(u,v,t, Δt)->norm(sys,u-v,Inf) # Time step error estimator
) where Tv
    inival=copy(inival)
    _initialize_dirichlet!(inival,system)
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
                    solve!(solution,inival, system ,control=control,tstep=Δt)
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
                        if !(control.force_first_step && istep==0)
                            throw(EmbeddingError())
                        else
                            solved=true
                        end
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
        sample(solution,times[i+1])
    end
    if control.verbose
        @printf("  Evolution: success\n")
    end
    return solution
end


"""
````
function solve(inival, system, times;
               control=NewtonControl(), 
               pre=function(sol,t) end,   
               post=function(sol,oldsol, t, Δt) end,
               sample=function(sol,t) end,
               delta=(u,v,t, Δt)->norm(sys,u-v,Inf),
               store_all=true,
               in_memory=true
)
````
Use implicit Euler method  + damped   Newton's method  to 
solve time dependent problem. Time step control is performed
according to the data in `control`.  All times in the vector
`times` are reached exactly.

Callbacks:
- `pre` is invoked before each time step
- `post`  is invoked after each time step
- `sample` is called for all times in `times[2:end]`.

`delta` is  used to control the time step.

If `store_all==true`, all timestep solutions are stored. Otherwise,
only solutions for the moments defined in `times` are stored.

Returns a transient solution object `sol` containing the stored solution,
see [`TransientSolution`](@ref)

"""
function solve(inival::AbstractMatrix,
               sys::AbstractSystem,
               times::AbstractVector;
               control=NewtonControl(),
               pre=function(sol,t) end,       # Function for preparing step
               post=function(sol,oldsol, t, Δt) end,      # Function for postprocessing successful step
               sample=function(sol,t) end,      # Function to be called for each t\in times[2:end]
               delta=(u,v,t, Δt)->norm(sys,u-v,Inf), # Time step error estimator
               store_all=true, # if true, store all solutions, otherwise, store only at sampling times
               in_memory=true
               )
    tsol=TransientSolution(Float64(times[1]),inival, in_memory=in_memory)
    solution=copy(inival)
    if store_all
        evolve!(solution,inival,sys,times,
                control=control,
                post=(u,uold,t,Δt)->(post(u,uold,t,Δt);append!(tsol,t,u)),
                pre=pre,
                sample=sample,
                delta=delta
                )
    else
        evolve!(solution,inival,sys,times,
                control=control,
                post=post,
                sample=(u,t)->(sample(u,t);append!(tsol,t,u)),
                pre=pre,
                delta=delta
                )
    end
    tsol
end


