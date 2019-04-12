

@inline function addnz(matrix,i,j,v::Tv,fac) where Tv
    if v!=zero(Tv)
        matrix[i,j]+=v*fac
    end
end

function _eval_and_assemble(this::AbstractSystem{Tv},
                            U::AbstractMatrix{Tv}, # Actual solution iteration
                            UOld::AbstractMatrix{Tv}, # Old timestep solution
                            F::AbstractMatrix{Tv},# Right hand side
                            tstep::Tv # time step size. Inf means stationary solution
                           ) where Tv

    if !isdefined(this,:matrix) # needed here for test function system 
        _create_matrix(this)
    end
    
    grid=this.grid

    physics::Physics=this.physics
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
        physics.flux(y,u,edge,physics.data)
    end

    function sourcewrap(y::AbstractVector)
        y.=0
        physics.source(y,node,physics.data)
    end

    function reactionwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        ## for ii in ..  uu[node.speclist[ii]]=u[ii]
        physics.reaction(y,u,node,physics.data)
        ## for ii in .. y[ii]=y[node.speclist[ii]]
    end

    function storagewrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.storage(y,u,node,physics.data)
    end

    @inline function breactionwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.breaction(y,u,bnode,physics.data)
    end

    @inline function bstoragewrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.bstorage(y,u,bnode,physics.data)
    end
    
    # Reset matrix + rhs
    matrix.nzval.=0.0
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
            fill!(node,grid, inode,icell)
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
                ForwardDiff.jacobian!(result_r,reactionwrap,Y,UK)
                res_react=DiffResults.value(result_r)
                jac_react=DiffResults.jacobian(result_r)
            end
            # Evaluate & differentiate storage term
            ForwardDiff.jacobian!(result_s,storagewrap,Y,UK)
            res_stor=DiffResults.value(result_s)
            jac_stor=DiffResults.jacobian(result_s)

            # Evaluate storage term for old timestep
            ### ALLOC ???
            storagewrap(oldstor,UKOld) 
            
            
            K=node.index

            for idof=firstnodedof(F,K):lastnodedof(F,K)
                ispec=spec(F,idof,K)
                add(F,idof,node_factors[inode]*(res_react[ispec]-src[ispec] + (res_stor[ispec]-oldstor[ispec])*tstepinv))
                for jdof=firstnodedof(F,K):lastnodedof(F,K)
                    jspec=spec(F,jdof,K)
                    addnz(matrix,idof,jdof,jac_react[ispec,jspec]+ jac_stor[ispec,jspec]*tstepinv,node_factors[inode])
                end
            end
        end
        
        for iedge=1:num_edges_per_cell(grid)
            if edge_factors[iedge]<edge_cutoff
                continue
            end

            fill!(edge,grid,iedge,icell)

            #Set up argument for fluxwrap
            for ispec=1:nspecies
                UKL[ispec]=U[ispec,edge.nodeK]
                UKL[ispec+nspecies]=U[ispec,edge.nodeL]
            end

            ForwardDiff.jacobian!(result_flx,fluxwrap,Y,UKL)
            res=DiffResults.value(result_flx)
            jac=DiffResults.jacobian(result_flx)

            K=edge.nodeK
            L=edge.nodeL
            fac=edge_factors[iedge]
            for idofK=firstnodedof(F,K):lastnodedof(F,K)
                ispec=spec(F,idofK,K)
                idofL=dof(F,ispec,L)
                if idofL==0
                    continue
                end
                add(F,idofK,fac*res[ispec])
                add(F,idofL,-fac*res[ispec])
                
                for jdofK=firstnodedof(F,K):lastnodedof(F,K)
                    jspec=spec(F,jdofK,K)
                    jdofL=dof(F,jspec,L)
                    if jdofL==0
                        continue
                    end
                    
                    addnz(matrix,idofK,jdofK,+jac[ispec,jspec            ],fac)
                    addnz(matrix,idofL,jdofK,-jac[ispec,jspec            ],fac)
                    addnz(matrix,idofK,jdofL,+jac[ispec,jspec+nspecies],fac)
                    addnz(matrix,idofL,jdofL,-jac[ispec,jspec+nspecies],fac)
                    
                end
            end
        end
    end

    for ibface=1:num_bfaces(grid)
        bfacefactors!(grid,ibface,bnode_factors)
        ibreg=grid.bfaceregions[ibface]

        for ibnode=1:num_nodes_per_bface(grid)
            fill!(bnode,grid,ibnode,ibface)
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
                        addnz(matrix,idof,idof,fac,1)
                    else
                        F[ispec,bnode.index]+=bnode_factors[ibnode]*(fac*U[ispec,bnode.index]-val)
                        addnz(matrix,idof,idof,fac,bnode_factors[ibnode])
                    end
                end
            end

            if isbreaction# involves bspecies and species
                ForwardDiff.jacobian!(result_r,breactionwrap,Y,UK)
                res_breact=DiffResults.value(result_r)
                jac_breact=DiffResults.jacobian(result_r)
                K=bnode.index
                fac=bnode_factors[ibnode]
                for idof=firstnodedof(F,K):lastnodedof(F,K)
                    ispec=spec(F,idof,K)
                    add(F,idof,fac*res_breact[ispec])
                    for jdof=firstnodedof(F,K):lastnodedof(F,K)
                        jspec=spec(F,jdof,K)
                        addnz(matrix,idof,jdof,jac_breact[ispec,jspec],fac)
                    end
                end

            end
            
            if isbstorage # should involve only bspecies
                # Evaluate & differentiate storage term
                ForwardDiff.jacobian!(result_s,bstoragewrap,Y,UK)
                res_bstor=DiffResults.value(result_s)
                jac_bstor=DiffResults.jacobian(result_s)
                
                # Evaluate storage term for old timestep
                bstoragewrap(oldbstor,UKOld)
                
                for idof=firstnodedof(F,K):lastnodedof(F,K)
                    ispec=spec(F,idof,K)
                    add(F,idof,fac*(res_bstor[ispec]-oldbstor[ispec])*tstepinv)
                    for jdof=firstnodedof(F,K):lastnodedof(F,K)
                        jspec=spec(F,jdof,K)
                        addnz(matrix,idof,jdof,jac_bstor[ispec,jspec],fac*tstepinv)
                    end
                end
            end
            
        end
    end
    _inactspecloop(this,U,UOld,F)
end



################################################################
function _solve!(
    solution::AbstractMatrix{Tv}, # old time step solution resp. initial value
    oldsol::AbstractMatrix{Tv}, # old time step solution resp. initial value
    this::AbstractSystem{Tv}, # Finite volume system
    control::NewtonControl,
    tstep::Tv
) where Tv
    
    if !isdefined(this,:matrix)
        _create_matrix(this)
        this.residual=unknowns(this)
        this.update=unknowns(this)
    end

    solution.=oldsol
    residual=this.residual
    update=this.update
    inidirichlet!(this,solution)

    # Newton iteration (quick and dirty...)
    oldnorm=1.0
    converged=false
    if control.verbose
        @printf("Start newton iteration: %s:%d\n", basename(@__FILE__),@__LINE__)
    end
    nlu=0
    lufact=nothing
    damp=control.damp_initial
    tolx=0.0

    for ii=1:control.max_iterations
        _eval_and_assemble(this,solution,oldsol,residual,tstep)
        
        # Sparse LU factorization
        # Here, we seem miss the possibility to re-use the 
        # previous symbolic information
        # We however reuse the factorization control.max_lureuse times.
        if nlu==0
            lufact=LinearAlgebra.lu(this.matrix)
            # LU triangular solve gives Newton update
            ldiv!(values(update),lufact,values(residual))
        else
            # When reusing lu factorization, we may try to iterate
            # Generally, this is advisable.
            if control.tol_linear <1.0
                bicgstabl!(values(update),this.matrix,values(residual),2,Pl=lufact,tol=control.tol_linear)
            else
                ldiv!(values(update),lufact,values(residual))
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
        if control.verbose
            @printf("  it=%03d norm=%.5e cont=%.5e\n",ii,norm, norm/oldnorm)
        end
        if norm<control.tol_absolute || norm <tolx
            converged=true
            break
        end
        oldnorm=norm
    end
    if !converged
        error("Error: no convergence")
    end
end

################################################################
"""
    function solve!(
        this::AbstractSystem{Tv},     # Finite volume system
        inival::AbstractMatrix{Tv},   # Initial value 
        solution::AbstractMatrix{Tv}; # Solution
        control=NewtonControl(),      # Newton solver control information
        tstep::Tv=Inf                 # Time step size. Inf means  stationary solution
    ) where Tv


Solution method for instance of System.

Perform solution of stationary system (if `tstep==Inf`) or implicit Euler time
step system. 

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
end



################################################################
"""
    function integrate(this::AbstractSystem{Tv},F::Function,U::AbstractMatrix{Tv}) where Tv

Integrate solution vector over domain. 
The results contains the integral for each species separately
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
            fill!(node,grid,inode,icell)
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

