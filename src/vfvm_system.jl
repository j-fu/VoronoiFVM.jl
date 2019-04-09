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

##########################################################
"""
   abstract type AbstractSystem

Abstract type for finite volume system structure
"""
abstract type AbstractSystem{Tv<:Number} end

abstract type Physics end

##################################################################
"""
    mutable struct SparseSystem{Tv}

Structure holding data for finite volume system solution.
Information on species distribution is kept in sparse
matrices, and the solution array is of type SparseSolutionArray,
i.e. effectively it is a sparse matrix.

Unlike in the DenseSystem, the system matrix handles exactly those
degrees of freedom which correspond to unknowns. However, handling
of the sparse matrix structures for the bookeeping of the unknowns
creates overhead.

"""
mutable struct SparseSystem{Tv} <: AbstractSystem{Tv}
    grid::Grid
    physics::Physics
    boundary_values::Array{Tv,2} # Array of boundary values  
    boundary_factors::Array{Tv,2}# Array of boundary factors 
    region_species::SparseMatrixCSC{Int8,Int16}
    bregion_species::SparseMatrixCSC{Int8,Int16}
    node_dof::SparseMatrixCSC{Int8,Int32}
    matrix::SparseMatrixCSC{Tv,Int32}
    species_homogeneous::Bool
    num_species::Int8
    SparseSystem{Tv}() where Tv = new()
end
##################################################################
"""
    function  SparseSystem(grid::Grid, physics::Physics, maxspec::Integer)

Constructor for SparseSystem. `physics` provides some user data, `maxspec`
is the maximum number of species.
"""
function  SparseSystem(grid::Grid,physics::Physics, maxspec::Integer)
    Tv=Base.eltype(grid)
    this=SparseSystem{Tv}()
    this.grid=grid
    this.physics=physics
    this.num_species=maxspec
    this.region_species=spzeros(Int8,Int16,maxspec,num_cellregions(grid))
    this.bregion_species=spzeros(Int8,Int16,maxspec,num_bfaceregions(grid))
    this.node_dof=spzeros(Int8,Int32,maxspec,num_nodes(grid))
    this.boundary_values=zeros(Tv,maxspec,num_bfaceregions(grid))
    this.boundary_factors=zeros(Tv,maxspec,num_bfaceregions(grid))
    this.species_homogeneous=false
    return this
end

##################################################################
"""
    mutable struct DenseSystem{Tv}

Structure holding data for finite volume system solution.
Information on species distribution is kept in dense
matrices, and the solution array is of type Array{2}.


Unlike in the SparseSystem, the system matrix handles exactly those
degrees of freedom which correspond to unknowns, and dummy 
degrees of freedom where unknowns are not defined. Handling
of the sparse matrix structures for the bookeeping of the unknowns
has less overhead, but additional dummy equations are added
to the system matrix.


"""
mutable struct DenseSystem{Tv} <: AbstractSystem{Tv}
    grid::Grid
    physics::Physics
    boundary_values::Array{Tv,2} # Array of boundary values  
    boundary_factors::Array{Tv,2}# Array of boundary factors 
    region_species::Array{Int8,2}
    bregion_species::Array{Int8,2}
    node_dof::Array{Int8,2}
    matrix::SparseMatrixCSC{Tv,Int32}
    species_homogeneous::Bool
    num_species::Int8
    DenseSystem{Tv}() where Tv = new()
end
##################################################################
"""
    function  DenseSystem(grid::Grid, physics::Physics, maxspec::Integer)

Constructor for DenseSystem. `physics` provides some user data, `maxspec`
is the maximum number of species.
"""
function  DenseSystem(grid::Grid,physics::Physics, maxspec::Integer)
    Tv=Base.eltype(grid)
    this=DenseSystem{Tv}()
    this.grid=grid
    this.physics=physics
    this.region_species=spzeros(Int8,Int16,maxspec,num_cellregions(grid))
    this.bregion_species=spzeros(Int8,Int16,maxspec,num_bfaceregions(grid))
    this.node_dof=spzeros(Int8,Int32,maxspec,num_nodes(grid))
    this.boundary_values=zeros(Tv,maxspec,num_bfaceregions(grid))
    this.boundary_factors=zeros(Tv,maxspec,num_bfaceregions(grid))
    this.species_homogeneous=false
    this.num_species=maxspec
    return this
end

##################################################################
"""
    function is_boundary_species(this::AbstractSystem, ispec::Integer)

Check if species number corresponds to boundary species.
"""
function is_boundary_species(this::AbstractSystem, ispec::Integer)
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
    function is_bulk_species(this::AbstractSystem, ispec::Integer)

Check if species number corresponds bulk species.
"""
function is_bulk_species(this::AbstractSystem, ispec::Integer)
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
    function add_species(this::AbstractSystem,ispec::Integer, regions::AbstractVector)

Add species to a list of bulk regions. Species numbers for
bulk and boundary species have to be distinct.
"""
function add_species(this::AbstractSystem,ispec::Integer, regions::AbstractVector)
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
    function add_boundary_species(this::AbstractSystem, ispec::Integer, regions::AbstractVector)

Add species to a list of boundary regions. Species numbers for
bulk and boundary species have to be distinct.

"""
function add_boundary_species(this::AbstractSystem, ispec::Integer, regions::AbstractVector)
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

# Create matrix in system
function _create_matrix(this::AbstractSystem)
    Tv=Base.eltype(this)
    this.matrix=spzeros(Tv,num_dof(this), num_dof(this))
    this.species_homogeneous=true
    for inode=1:size(this.node_dof,2)
        for ispec=1:size(this.node_dof,1)
            if this.node_dof[ispec,inode]!=ispec
                this.species_homogeneous=false
                return
            end
        end
    end
end

##################################################################
"""
    num_dof(this::SparseSystem)

Number of degrees of freedom for system.
"""
num_dof(this::SparseSystem)= nnz(this.node_dof)

##################################################################
"""
    num_dof(this::SparseSystem)

Number of degrees of freedom for system.
"""
num_dof(this::DenseSystem)= length(this.node_dof)

##################################################################
"""
    num_species(this::SparseSystem)

Number of species in system
"""
num_species(this::AbstractSystem{Tv}) where Tv = this.num_species





##################################################################
"""
    struct SparseSolutionArray{Tv} <: AbstractMatrix{Tv}
        node_dof::SparseMatrixCSC{Tv,Int16}
    end

Struct holding solution information for SparseSystem. Solution
is stored in a sparse matrix structure.

This class plays well with the abstract array interface.
"""
struct SparseSolutionArray{Tv} <: AbstractMatrix{Tv}
    node_dof::SparseMatrixCSC{Tv,Int16}
end

##################################################################
"""
    function unknowns(::SparseSystem)

Create a solution vector for system.
"""
function unknowns(sys::SparseSystem{Tv}) where Tv
    return SparseSolutionArray{Tv}(SparseMatrixCSC(sys.node_dof.m,
                                        sys.node_dof.n,
                                        sys.node_dof.colptr,
                                        sys.node_dof.rowval,
                                        Array{Tv}(undef,num_dof(sys))
                                        )
                        )
end


##################################################################
"""
    function unknowns(::DenseSystem)

Create a solution vector for system.
"""
unknowns(sys::DenseSystem{Tv}) where Tv=Array{Tv}(undef,size(sys.node_dof,1), size(sys.node_dof,2))


##################################################################
"""
    size(a::SparseSolutionArray)
    
Return size of solution array.
"""
Base.size(a::SparseSolutionArray)=size(a.node_dof)

##################################################################
"""
    num_nodes(a)
                        
Number of nodes (size of second dimension) of solution array.
"""
num_nodes(a)=size(a,2)

##################################################################
"""
    num_species(a)

Number of species (size of first dimension) of solution array.
"""
num_species(a)=size(a,1)

##################################################################
"""
    values(a::SparseSolutionArray)

Array of values in solution array.
"""
values(a::SparseSolutionArray)=a.node_dof.nzval


##################################################################
"""
    values(a::Array)

Array of values in solution array.
"""
values(a::Array)= vec(a)



##################################################################
"""
    copy(this::SparseSolutionArray)

Create a copy of solution array
"""
Base.copy(this::SparseSolutionArray{Tv}) where Tv = SparseSolutionArray{Tv}(SparseMatrixCSC(this.node_dof.m,
                                                                                            this.node_dof.n,
                                                                                            this.node_dof.colptr,
                                                                                            this.node_dof.rowval,
                                                                                            Base.copy(this.node_dof.nzval)
                                                                                            )
                                                                            )
##################################################################
"""
    function dof(a::SparseSolutionArray,ispec, inode)

Get number of degree of freedom. Return 0 if species is not defined in node.
"""
@inline function dof(a::SparseSolutionArray{Tv},i::Integer, j::Integer) where Tv
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
    function dof(Array{2},ispec, inode)

Get number of degree of freedom.
"""
dof(a::Array{Tv,2}, ispec::Integer, K::Integer) where Tv = (K-1)*size(a,1)+ispec


##################################################################
"""
    function setdof!(a::SparseSolutionArray,v,i::Integer)

Set value for degree of freedom.
"""
function setdof!(a::SparseSolutionArray,v,i::Integer)
    a.node_dof.nzval[i] = v
end

##################################################################
"""
    function getdof(a::SparseSolutionArray,i::Integer)

Return  value for degree of freedom.
"""
getdof(a::SparseSolutionArray,i::Integer) =a.node_dof.nzval[i] 




##################################################################
"""
     setindex!(a::SparseSolutionArray, v, ispec, inode)

Accessor for solution array.
"""
function Base.setindex!(a::SparseSolutionArray, v, ispec::Integer, inode::Integer)
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
     getindex!(a::SparseSolutionArray, ispec, inode)

Accessor for solution array.
"""
function Base.getindex(a::SparseSolutionArray, ispec::Integer, inode::Integer)
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
    struct SubgridArrayView{Tv} <: AbstractMatrix{Tv}

Struct holding information for solution array view on subgrid
"""
struct SubgridArrayView{Tv,Ta} <: AbstractMatrix{Tv}
    sysarray::Ta
    subgrid::SubGrid
end

##################################################################
"""
    view(a::AbstractMatrix{Tv},sg::SubGrid)

Create a view of the solution array on a subgrid.
"""
Base.view(a::AbstractMatrix{Tv},sg::SubGrid) where Tv = SubgridArrayView{Tv,typeof(a)}(a,sg)


##############################################################################
"""
    getindex(aview::SubgridArrayView,ispec::Integer,inode::Integer)

Accessor method for subgrid array view.
"""
Base.getindex(aview::SubgridArrayView,ispec::Integer,inode::Integer) = aview.sysarray[ispec,aview.subgrid.node_in_parent[inode]]

##############################################################################
"""
    setindex!(aview::SubgridArrayView,v,ispec::Integer,inode::Integer)

Accessor method for subgrid array view.
"""
@inline function Base.setindex!(aview::SubgridArrayView,v,ispec::Integer,inode::Integer)
    aview.sysarray[ispec,aview.subgrid.node_in_parent[inode]]=v
    return aview
end

##################################################################
"""
    size(a::SubgridArrayView)
    
Return size of solution array view.
"""
Base.size(a::SubgridArrayView)=(size(a.sysarray,1),size(a.subgrid.node_in_parent,1))



isdof(this::AbstractSystem,ispec,inode)= this.node_dof[ispec,inode]==ispec ? true : false






##############################################################################
"""
    function inidirichlet!(this::AbstractSystem,U)

Initialize dirichlet boundary values for solution.
"""
function inidirichlet!(this::AbstractSystem,U) where Tv
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
    _inactspecinit(this,U)
end


#############################################################################
# Assemble routine

function _inactspecloop(this::DenseSystem,U,Uold,F)
    if this.species_homogeneous
        return
    end
    for inode=1:size(this.node_dof,2)
        for ispec=1:size(this.node_dof,1)
            if !isdof(this,ispec,inode)
                F[ispec,inode]+= U[ispec,inode]-Uold[ispec,inode];
                idof=dof(F,ispec,inode)
                this.matrix[idof,idof]+=1.0
            end
        end
    end
end

function _inactspecinit(this::DenseSystem,U)
    if this.species_homogeneous
        return
    end
    for inode=1:size(this.node_dof,2)
        for ispec=1:size(this.node_dof,1)
            if !isdof(this,ispec,inode)
                U[ispec,inode]=0
            end
        end
    end
end


function _inactspecloop(this::SparseSystem,U,Uold,F)
end

function _inactspecinit(this::SparseSystem,U)
end


firstnodedof(U::SparseSolutionArray{Tv},K) where Tv =U.node_dof.colptr[K]
lastnodedof(U::SparseSolutionArray{Tv},K) where Tv=U.node_dof.colptr[K+1]-1
spec(U::SparseSolutionArray{Tv},idof,K) where Tv=U.node_dof.rowval[idof]
add(U::SparseSolutionArray{Tv},idof,val) where Tv=U.node_dof.nzval[idof]+=val

firstnodedof(U::Matrix{Tv},K) where Tv = (K-1)*size(U,1)+1
lastnodedof(U::Matrix{Tv},K) where Tv = K*size(U,1)
spec(U::Matrix{Tv},idof,K) where Tv =   idof-(K-1)*size(U,1)
add(U::Matrix{Tv},idof,val) where Tv=vec(U)[idof]+=val




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

    grid=this.grid

    physics::Physics=this.physics
    node::Node=Node{Tv}()
    bnode::BNode=BNode{Tv}()
    edge::Edge=Edge{Tv}()
    edge_cutoff=1.0e-12
    nspecies::Int32=num_species(this)
    
    if !isdefined(this,:matrix)
        _create_matrix(this)
    end
    matrix=this.matrix
    
    
    K1::Int32=1
    KN::Int32=nspecies
    L1::Int32=nspecies+1
    LN::Int32=2*nspecies

    @inline function fluxwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        @views physics.flux(physics,edge,y,u[K1:KN],u[L1:LN])
    end
    
    @inline function sourcewrap(y::AbstractVector)
        y.=0
        physics.source(physics,node,y)
    end

    @inline function reactionwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        ## for ii in ..  uu[node.speclist[ii]]=u[ii]
        physics.reaction(physics,node,y,u)
        ## for ii in .. y[ii]=y[node.speclist[ii]]
    end

    @inline function storagewrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.storage(physics,node,y,u)
    end

    @inline function breactionwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.breaction(physics,bnode,y,u)
    end

    @inline function bstoragewrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.bstorage(physics,bnode,y,u)
    end
    
    # Reset matrix + rhs
    matrix.nzval.=0.0
    F.=0.0

    # structs holding diff results for storage, reaction,  flux ...
    result_r=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_s=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_br=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_bs=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
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
            @views begin
                # xx gather:
                # ii=0
                # for i region_spec.colptr[ireg]:region_spec.colptr[ireg+1]-1
                #    ispec=Fdof.rowval[idof]
                #    ii=ii+1
                #    node.speclist[ii]=ispec
                #    UK[ii]=U[ispec,K]
                UK[1:nspecies]=U[:,node.index]
                UKOld[1:nspecies]=UOld[:,node.index]
                # Evaluate source term
            end

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

            @views begin
                fill!(edge,grid,iedge,icell)

                #Set up argument for fluxwrap
                UKL[K1:KN]=U[:,edge.nodeK]
                UKL[L1:LN]=U[:,edge.nodeL]
                
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
                    addnz(matrix,idofK,jdofL,+jac[ispec,jspec+nspecies],fac)
                    addnz(matrix,idofL,jdofK,-jac[ispec,jspec            ],fac)
                    addnz(matrix,idofL,jdofL,-jac[ispec,jspec+nspecies],fac)
                    
                end
            end
            
        end
    end

    for ibface=1:num_bfaces(grid)
        bfacefactors!(grid,ibface,bnode_factors)
        ibreg=grid.bfaceregions[ibface]
        bnode.region=ibreg
        for ibnode=1:num_nodes_per_bface(grid)
            @views begin
                fill!(bnode,grid,ibnode,ibface)

                UK[1:nspecies]=U[:,bnode.index]
                UKOld[1:nspecies]=UOld[:,bnode.index]
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

            if isdefined(physics, :breaction)# involves bspecies and species
                ForwardDiff.jacobian!(result_br,breactionwrap,Y,UK)
                res_breact=DiffResults.value(result_br)
                jac_breact=DiffResults.jacobian(result_br)
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
            
            if isdefined(physics, :bstorage) # should involve only bspecies
                # Evaluate & differentiate storage term
                ForwardDiff.jacobian!(result_bs,bstoragewrap,Y,UK)
                res_bstor=DiffResults.value(result_bs)
                jac_bstor=DiffResults.jacobian(result_bs)
                
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
function _solve(
    this::AbstractSystem{Tv}, # Finite volume system
    oldsol::AbstractMatrix{Tv}, # old time step solution resp. initial value
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
        oldsol::AbstractMatrix;     # old time step solution resp. initial value
        control=NewtonControl(), # Solver control information (optional)
        tstep::Tv=Inf            # Time step size. Inf means  stationary solution. (optional)
        )


Solution method for instance of System.

Perform solution of stationary system (if `tstep==Inf`) or implicit Euler time
step system. 

"""
function solve(
    this::AbstractSystem{Tv}, # Finite volume system
    oldsol::AbstractMatrix{Tv}; # old time step solution resp. initial value
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
function integrate(this::AbstractSystem,F::Function,U)
````

Integrate solution vector over domain. Returns an `Array{Tv,1}`
containing the integral for each species.
"""
function integrate(this::AbstractSystem{Tv},F::Function,U::AbstractMatrix{Tv}) where Tv
    grid=this.grid
    nspecies=num_species(this)
    integral=zeros(Tv,nspecies)
    res=zeros(Tv,nspecies)
    node=Node{Tv}()
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))

    for icell=1:num_cells(grid)
        cellfactors!(grid,icell,node_factors,edge_factors)
        for inode=1:num_nodes_per_cell(grid)
            fill!(node,grid,inode,icell)
            F(this.physics,node,res,U[:,node.index])
            for ispec=1:nspecies
                if this.node_dof[ispec,node.index]==ispec
                    integral[ispec]+=node_factors[inode]*res[ispec]
                end
            end
        end
    end
    return integral
end

