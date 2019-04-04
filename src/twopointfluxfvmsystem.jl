# Packages for Autodiff magic. These need to be installed via Pkg
using ForwardDiff, DiffResults
using IterativeSolvers

# These are in the standard distro
using SparseArrays
using LinearAlgebra
using Printf

##################################################################
"""
Constant to be used as boundary condition factor 
to mark Dirichlet boundary conditons.    
"""
const Dirichlet=1.0e30

##################################################################
"""
    value(x)

Extract value from dual number. Use to debug physics callbacks.
Re-exported from ForwardDiff.jl
"""
const value=ForwardDiff.value

##################################################################
"""
    mutable struct System{Tv}

Main structure holding data for system solution.

"""
mutable struct System{Tv}
    grid::Grid
    physics::Any
    boundary_values::Array{Tv,2} # Array of boundary values  
    boundary_factors::Array{Tv,2}# Array of boundary factors 
    region_species::SparseMatrixCSC{Int8,Int16}
    bregion_species::SparseMatrixCSC{Int8,Int16}
    node_dof::SparseMatrixCSC{Int8,Int32}
    matrix::SparseArrays.SparseMatrixCSC{Tv,Int32}
    function System{Tv}() where Tv
        return new{Tv}()
    end
end
##################################################################
"""
    function  System(grid::Grid, physics::Any, maxspec::Integer)

Constructor for System. `physics` provides some user data, `maxspec`
is the maximum number of species.
"""
function  System(grid::Grid,physics, maxspec::Integer)
    Tv=Base.eltype(grid)
    this=System{Tv}()
    this.grid=grid
    this.physics=physics
    this.region_species=spzeros(Int8,Int16,maxspec,num_cellregions(grid))
    this.bregion_species=spzeros(Int8,Int16,maxspec,num_bfaceregions(grid))
    this.node_dof=spzeros(Int8,Int32,maxspec,num_nodes(grid))
    this.boundary_values=zeros(Tv,maxspec,num_bfaceregions(grid))
    this.boundary_factors=zeros(Tv,maxspec,num_bfaceregions(grid))
    return this
end

##################################################################
"""
    function is_boundary_species(this::System, ispec::Integer)

Check if species number corresponds to boundary species.
"""
function is_boundary_species(this::System, ispec::Integer)
    isbspec=false
    for ibreg=1:num_bfaceregions(this.grid)
        if this.bregion_species[ispec,ibreg]>0
            isbspec=true
        end
    end
    return isbspec
end

##################################################################
"""
    function is_bulk_species(this::System, ispec::Integer)

Check if species number corresponds bulk species.
"""
function is_bulk_species(this::System, ispec::Integer)
    isrspec=false
    for ixreg=1:num_cellregions(this.grid)
        if this.region_species[ispec,ixreg]>0
            isrspec=true
        end
    end
    return isrspec
end

##################################################################
"""
    function add_species(this::System,ispec::Integer, regions::AbstractArray)

Add species to a list of bulk regions. Species numbers for
bulk and boundary species have to be distinct.
"""
function add_species(this::System,ispec::Integer, regions::AbstractArray)
    if is_boundary_species(this,ispec)
        throw(DomainError(ispec,"Species is already boundary species"))
    end

    for i in eachindex(regions)
        ireg=regions[i]
        this.region_species[ispec,ireg]=ispec
        for icell=1:num_cells(this.grid)
            if this.grid.cellregions[icell]==ireg
                for inode=1:size(this.grid.cellnodes,1)
                    this.node_dof[ispec,this.grid.cellnodes[inode,icell]]=ispec
                end
            end
        end
    end
end

##################################################################
"""
    function add_boundary_species(this::System, ispec::Integer, regions::AbstractArray)

Add species to a list of boundary regions. Species numbers for
bulk and boundary species have to be distinct.

"""
function add_boundary_species(this::System, ispec::Integer, regions::AbstractArray)
    if is_bulk_species(this,ispec)
        throw(DomainError(ispec,"Species is already bulk species"))
    end
    for i in eachindex(regions)
        ireg=regions[i]
        this.bregion_species[ispec,ireg]=1
        for ibface=1:num_bfaces(this.grid)
            if this.grid.bfaceregions[ibface]==ireg
                for inode=1:size(this.grid.bfacenodes,1)
                    this.node_dof[ispec,this.grid.bfacenodes[inode,ibface]]=ispec
                end
            end
        end
    end
end

##################################################################
"""
    num_dof(this::System)

Number of degrees of freedom for system.
"""
num_dof(this::System)= nnz(this.node_dof)

##################################################################
"""
    num_species(this::System)

Number of species in system
"""
num_species(this::System)= this.node_dof.m








##################################################################
"""
    struct SysArray{Tv} <: AbstractArray{Tv,2}
        node_dof::SparseMatrixCSC{Tv,Int16}
    end

Struct holding solution information for system. Solution
is stored in a sparse matrix structure.

This class plays well with the abstract array interface
"""
struct SysArray{Tv} <: AbstractArray{Tv,2}
    node_dof::SparseMatrixCSC{Tv,Int16}
end

##################################################################
"""
    function unknowns(system)

Create a solution vector for system.
"""
function unknowns(sys::System{Tv}) where Tv
    return SysArray{Tv}(SparseMatrixCSC(sys.node_dof.m,
                                        sys.node_dof.n,
                                        sys.node_dof.colptr,
                                        sys.node_dof.rowval,
                                        Array{Tv}(undef,num_dof(sys))
                                        )
                        )
end

##################################################################
"""
    size(a::SysArray)
    
Return size of solution array.
"""
Base.size(a::SysArray)=size(a.node_dof)

##################################################################
"""
    num_nodes(a::SysArray)
                        
Number of nodes (size of second dimension) of solution array.
"""
num_nodes(a::SysArray)=size(a,2)

##################################################################
"""
    num_species(a::SysArray)

Number of species (size of first dimension) of solution array.
"""
num_species(a::SysArray)=size(a,1)

##################################################################
"""
    values(a::SysArray)=a.node_dof

Array of values in solution array.
"""
values(a::SysArray)=a.node_dof.nzval


##################################################################
"""
    copy(this::SysArray)

Create a copy of solution array
"""
function Base.copy(this::SysArray{Tv}) where Tv
    return SysArray{Tv}(SparseMatrixCSC(this.node_dof.m,
                                        this.node_dof.n,
                                        this.node_dof.colptr,
                                        this.node_dof.rowval,
                                        Base.copy(this.node_dof.nzval)
                                        )
                        )
end

##################################################################
"""
    function dof(a::SysArray,ispec, inode)

Get number of degree of freedom. Return 0 if species is not defined in node.
"""
@inline function dof(a::SysArray{Tv},i::Integer, j::Integer) where Tv
    A=a.node_dof
    coljfirstk = Int(A.colptr[j])
    coljlastk = Int(A.colptr[j+1] - 1)
    searchk = searchsortedfirst(A.rowval, i, coljfirstk, coljlastk, Base.Order.Forward)
    if searchk <= coljlastk && A.rowval[searchk] == i
        return searchk
    end
    return 0
end

##################################################################
"""
    function setdof!(a::SysArray,v,i::Integer)

Set value for degree of freedom.
"""
@inline function setdof!(a::SysArray,v,i::Integer)
    a.node_dof.nzval[i] = v
end

##################################################################
"""
    function getdof(a::SysArray,i::Integer)

Return  value for degree of freedom.
"""
@inline function getdof(a::SysArray,i::Integer)
    return a.node_dof.nzval[i] 
end

##################################################################
"""
     setindex!(a::SysArray, v, ispec, inode)

Accessor for solution array.
"""
function Base.setindex!(a::SysArray, v, ispec::Integer, inode::Integer)
    searchk=dof(a,ispec,inode)
    if searchk>0
        setdof!(a,v,searchk)
        return a
    end
    # TODO: what is the right reacton here ?
    # Ignoring seems to be better, so we can broacast etc.
    # throw(DomainError("undefined degree of freedom"))
end

##################################################################
"""
     getindex!(a::SysArray, ispec, inode)

Accessor for solution array.
"""
function Base.getindex(a::SysArray, ispec::Integer, inode::Integer)
    searchk=dof(a,ispec,inode)
    if searchk>0
        return getdof(a,searchk)
    end
    #
    # TODO: what is the right reacton here ?
    # Actually, NaN plays well with pyplot...
    return NaN
end

##################################################################
"""
    struct SubgridSysArrayView{Tv} <: AbstractArray{Tv,2}

Struct holding information for solution array view on subgrid
"""
struct SubgridSysArrayView{Tv} <: AbstractArray{Tv,2}
    sysarray::SysArray{Tv}
    subgrid::SubGrid
end

##################################################################
"""
    view(a::SysArray{Tv},sg::SubGrid)

Create a view of the solution array on a subgrid.
"""
function Base.view(a::SysArray{Tv},sg::SubGrid) where Tv
    return SubgridSysArrayView{Tv}(a,sg)
end

##############################################################################
"""
    getindex(aview::SubgridSysArrayView,ispec::Integer,inode::Integer)

Accessor method for subgrid array view.
"""
function Base.getindex(aview::SubgridSysArrayView,ispec::Integer,inode::Integer)
    return aview.sysarray[ispec,aview.subgrid.node_in_parent[inode]]
end

##############################################################################
"""
    setindex!(aview::SubgridSysArrayView,v,ispec::Integer,inode::Integer)

Accessor method for subgrid array view.
"""
function Base.setindex!(aview::SubgridSysArrayView,v,ispec::Integer,inode::Integer)
    aview.sysarray[ispec,aview.subgrid.node_in_parent[inode]]=v
    return aview
end

##################################################################
"""
    size(a::SubgridSysArrayView)
    
Return size of solution array view.
"""
Base.size(a::SubgridSysArrayView)=(size(a.sysarray,1),size(a.subgrid.node_in_parent,1))







##############################################################################
"""
    function inidirichlet!(this::System,U)

Initialize dirichlet boundary values for solution.
"""
function inidirichlet!(this::System{Tv},U::SysArray{Tv}) where Tv
    for ibface=1:num_bfaces(this.grid)
        ibreg=this.grid.bfaceregions[ibface]
        for ispec=1:num_species(this)
            if this.boundary_factors[ispec,ibreg]==Dirichlet
                for inode=1:dim_grid(this.grid)
                    U[ispec,this.grid.bfacenodes[inode,ibface]]=this.boundary_values[ispec,ibreg]
                end
            end
        end
    end
end


#############################################################################
# Assemble routine
function _eval_and_assemble(this::System,
                            U, # Actual solution iteration
                            UOld, # Old timestep solution
                            F,# Right hand side
                            tstep # time step size. Inf means stationary solution
                           )

    grid=this.grid
    Tv=Base.eltype(grid)

    physics=this.physics
    node=Node()
    edge=Edge()
    edge_cutoff=1.0e-12
    nspecies=num_species(this)
    
    
    if !isdefined(this,:matrix)
        this.matrix=spzeros(Tv,num_dof(this), num_dof(this))
    end

    @inline function addnz(matrix,i,j,v,fac)
        if v!=0.0
            matrix[i,j]+=v*fac
        end
    end
    
    
    K1=1
    KN=nspecies
    L1=nspecies+1
    LN=2*nspecies

    @inline function fluxwrap(y,u)
        y.=0
        @views physics.flux(physics,node,y,u[K1:KN],u[L1:LN])
    end
    
    @inline function sourcewrap(y)
        y.=0
        physics.source(physics,node,y)
    end

    @inline function reactionwrap(y,u)
        y.=0
        physics.reaction(physics,node,y,u)
    end

    @inline function storagewrap(y,u)
        y.=0
        physics.storage(physics,node,y,u)
    end

    @inline function breactionwrap(y,u)
        y.=0
        physics.breaction(physics,node,y,u)
    end

    @inline function bstoragewrap(y,u)
        y.=0
        physics.bstorage(physics,node,y,u)
    end
    
    M=this.matrix
    
    # Reset matrix + rhs
    M.nzval.=0.0
    F.=0.0

    # structs holding diff results for storage, reaction,  flux ...
    result_r=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_s=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_br=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_bs=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_flx=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,2*nspecies))

    # Arrays holding function results
    Y=Array{Tv,1}(undef,nspecies)
    res_react=zeros(Tv,nspecies)
    jac_react=zeros(Tv,nspecies,nspecies)

    # Arrays for gathering solution data
    UK=Array{Tv,1}(undef,nspecies)
    UKOld=Array{Tv,1}(undef,nspecies)
    UKL=Array{Tv,1}(undef,2*nspecies)

    # array holding source term
    src=zeros(Tv,nspecies)

    # arrays holding storage terms for old solution
    oldstor=zeros(Tv,nspecies)
    oldbstor=zeros(Tv,nspecies)


    # Inverse of timestep
    # According to Julia documentation, 1/Inf=0 which
    # comes handy to write compact code here.
    tstepinv=1.0/tstep 

    
    # Arrays holding for factor data
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))

    # Main cell loop
    for icell=1:num_cells(grid)
        # set up form factors
        cellfactors!(grid,icell,node_factors,edge_factors)

        # set up data for callbacks
        node.region=reg_cell(grid,icell)
        edge.region=reg_cell(grid,icell)

        
        for inode=1:num_nodes_per_cell(grid)
            @views begin
                K=cellnode(grid,inode,icell)
                node.index=K
                node.coord=nodecoord(grid,K)
                UK[1:nspecies]=U[:,K]
                UKOld[1:nspecies]=UOld[:,K]
            end
            # Evaluate source term
            if isdefined(physics,:source)
                sourcewrap(src)
            end
            
            # Evaluate & differentiate storage term
            ForwardDiff.jacobian!(result_s,storagewrap,Y,UK)
            res_stor=DiffResults.value(result_s)
            jac_stor=DiffResults.jacobian(result_s)
            
            # Evaluate storage term for old timestep
            storagewrap(oldstor,UKOld)
            
            # Evaluate reaction term if present
            if isdefined(physics, :reaction)
                ForwardDiff.jacobian!(result_r,reactionwrap,Y,UK)
                res_react=DiffResults.value(result_r)
                jac_react=DiffResults.jacobian(result_r)
            end
            
            # Assembly results and jacobians
            Fdof=F.node_dof
            for idof=Fdof.colptr[K]:Fdof.colptr[K+1]-1
                ispec=Fdof.rowval[idof]
                Fdof.nzval[idof]+=node_factors[inode]*(res_react[ispec]-src[ispec] + (res_stor[ispec]-oldstor[ispec])*tstepinv)
                for jdof=Fdof.colptr[K]:Fdof.colptr[K+1]-1
                    jspec=Fdof.rowval[jdof]
                    addnz(M,idof,jdof,jac_react[ispec,jspec]+ jac_stor[ispec,jspec]*tstepinv,node_factors[inode])
                end
            end
        end
        
        for iedge=1:num_edges_per_cell(grid)
            if edge_factors[iedge]<edge_cutoff
                continue
            end
            @views begin
                K=celledgenode(grid,1,iedge,icell)
                L=celledgenode(grid,2,iedge,icell)
                edge.index=iedge
                edge.nodeK=K
                edge.nodeL=L
                edge.coordL=nodecoord(grid,L)
                edge.coordK=nodecoord(grid,K)
                
                #Set up argument for fluxwrap
                UKL[K1:KN]=U[:,K]
                UKL[L1:LN]=U[:,L]
            end
                
            ForwardDiff.jacobian!(result_flx,fluxwrap,Y,UKL)
            
            res=DiffResults.value(result_flx)
            jac=DiffResults.jacobian(result_flx)


            Fdof=F.node_dof
            for idofK=Fdof.colptr[K]:Fdof.colptr[K+1]-1
                ispec=Fdof.rowval[idofK]
                idofL=dof(F,ispec,L)
                if idofL==0
                    continue
                end

                Fdof.nzval[idofK]+=res[ispec]*edge_factors[iedge]
                Fdof.nzval[idofL]-=res[ispec]*edge_factors[iedge]

                for jdofK=Fdof.colptr[K]:Fdof.colptr[K+1]-1
                    jspec=Fdof.rowval[jdofK]
                    jdofL=dof(F,jspec,L)
                    if jdofL==0
                        continue
                    end
                    
                    addnz(M,idofK,jdofK,+jac[ispec,jspec            ],edge_factors[iedge])
                    addnz(M,idofK,jdofL,+jac[ispec,jspec+nspecies],edge_factors[iedge])
                    addnz(M,idofL,jdofK,-jac[ispec,jspec            ],edge_factors[iedge])
                    addnz(M,idofL,jdofL,-jac[ispec,jspec+nspecies],edge_factors[iedge])
                    
                end
            end
        end
    end

   bnode_factors=zeros(Tv,num_nodes_per_bface(grid))
   for ibface=1:num_bfaces(grid)
        bfacefactors!(grid,ibface,bnode_factors)
        ibreg=grid.bfaceregions[ibface]
        node.region=ibreg
        for ibnode=1:num_nodes_per_bface(grid)
            @views begin
                K=bfacenode(grid,ibnode,ibface)
                node.index=K
                node.coord=nodecoord(grid,K)
                UK[1:nspecies]=U[:,K]
                UKOld[1:nspecies]=UOld[:,K]
            end

            for ispec=1:nspecies # should involve only rspecies
                fac=this.boundary_factors[ispec,ibreg]
                val=this.boundary_values[ispec,ibreg]
                if fac!=Dirichlet
                    fac*=bnode_factors[ibnode]
                end
                idof=dof(F,ispec,K)
                if idof>0
                    F[ispec,K]+=fac*(U[ispec,K]-val)
                    addnz(M,idof,idof,fac,1)
                end
            end
            
            if isdefined(physics, :breaction)# involves bspecies and species
                ForwardDiff.jacobian!(result_br,breactionwrap,Y,UK)
                res_breact=DiffResults.value(result_br)
                jac_breact=DiffResults.jacobian(result_br)
                Fdof=F.node_dof
                for idof=Fdof.colptr[K]:Fdof.colptr[K+1]-1
                    ispec=Fdof.rowval[idof]
                    Fdof.nzval[idof]+=bnode_factors[ibnode]*res_breact[ispec]
                    for jdof=Fdof.colptr[K]:Fdof.colptr[K+1]-1
                        jspec=Fdof.rowval[jdof]
                        addnz(M,idof,jdof, jac_breact[ispec,jspec],bnode_factors[ibnode])
                    end
                end
            end
            
            if isdefined(physics, :bstorage) # should involve only bspecies
                # Evaluate & differentiate storage term
                ForwardDiff.jacobian!(result_bs,bstoragewrap,Y,UK)
                res_bstor=DiffResults.value(result_bs)
                jac_bstor=DiffResults.jacobian(result_bs)
                
                # Evaluate storage term for old timestep
                bstoragewrap(oldbstor,UKOld)
                
                Fdof=F.node_dof
                for idof=Fdof.colptr[K]:Fdof.colptr[K+1]-1
                    ispec=Fdof.rowval[idof]
                    Fdof.nzval[idof]+=bnode_factors[ibnode]*(res_bstor[ispec]-oldbstor[ispec])*tstepinv
                    for jdof=Fdof.colptr[K]:Fdof.colptr[K+1]-1
                        jspec=Fdof.rowval[jdof]
                        if jac_bstor[ispec,jspec]==0.0
                            continue
                        end
                        addnz(M,idof,jdof,jac_bstor[ispec,jspec],bnode_factors[ibnode]*tstepinv)
                    end
                end
            end
        end
   end

end



################################################################
function _solve(
    this::System{Tv}, # Finite volume system
    oldsol::SysArray{Tv}, # old time step solution resp. initial value
    control::NewtonControl,
    tstep::Tv
) where Tv
    
    solution=copy(oldsol)
    residual=copy(solution)
    update=copy(solution)
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
    return solution
end

################################################################
"""
    function solve(
        this::System,            # Finite volume system
        oldsol::Array{Tv,1};     # old time step solution resp. initial value
        control=NewtonControl(), # Solver control information (optional)
        tstep::Tv=Inf            # Time step size. Inf means  stationary solution. (optional)
        )


Solution method for instance of System.

Perform solution of stationary system (if `tstep==Inf`) or implicit Euler time
step system. 

"""
function solve(
    this::System{Tv}, # Finite volume system
    oldsol::SysArray{Tv}; # old time step solution resp. initial value
    control=NewtonControl(), # Newton solver control information
    tstep::Tv=Inf          # Time step size. Inf means  stationary solution
) where Tv
    if control.verbose
        @time begin
            retval= _solve(this,oldsol,control,tstep)
        end
        return retval
    else
        return _solve(this,oldsol,control,tstep)
    end
end



"""
````
function integrate(this::System,F::Function,U)
````

Integrate solution vector over domain. Returns an `Array{Int64,1}`
containing the integral for each species.
"""
function integrate(this::System{Tv},F::Function,U::SysArray{Tv}) where Tv
    grid=this.grid
    nspecies=num_species(this)
    integral=zeros(Tv,nspecies)
    res=zeros(Tv,nspecies)
    node=Node()
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))

    for icell=1:num_cells(grid)
        cellfactors!(grid,icell,node_factors,edge_factors)
        node.region=reg_cell(grid,icell)

        for inode=1:num_nodes_per_cell(grid)
            K=cellnode(grid,inode,icell)
            node.index=K
            node.coord=nodecoord(grid,K)
            F(this.physics,node,res,U[:,K])
            for ispec=1:nspecies
                integral[ispec]+=node_factors[inode]*res[ispec]
            end
        end
    end
    return integral
end

