##########################################################
"""
       abstract type AbstractSystem
    
Abstract type for finite volume system structure
"""
abstract type AbstractSystem{Tv<:Number} end


##################################################################
"""
Constant to be used as boundary condition factor 
to mark Dirichlet boundary conditons.    
"""
const Dirichlet=1.0e30



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
    update::SparseSolutionArray{Tv}
    residual::SparseSolutionArray{Tv}
    SparseSystem{Tv}() where Tv = new()
end
##################################################################
"""
    function  SparseSystem(grid::Grid, physics::Physics, maxspec::Integer)

Constructor for SparseSystem. `physics` provides some user data, `maxspec`
is the maximum number of species.
"""
function  SparseSystem(grid::Grid,physics::Physics)
    Tv=Base.eltype(grid)
    this=SparseSystem{Tv}()
    maxspec=physics.num_species
    this.grid=grid
    this.physics=physics
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
    update::Matrix{Tv}
    residual::Matrix{Tv}
    DenseSystem{Tv}() where Tv = new()
end
##################################################################
"""
    function  DenseSystem(grid::Grid, physics::Physics, maxspec::Integer)

Constructor for DenseSystem. `physics` provides some user data, `maxspec`
is the maximum number of species.
"""
function  DenseSystem(grid::Grid,physics::Physics)
    Tv=Base.eltype(grid)
    this=DenseSystem{Tv}()
    maxspec=physics.num_species
    this.grid=grid
    this.physics=physics
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
    function enable_species(this::AbstractSystem,ispec::Integer, regions::AbstractVector)

Add species to a list of bulk regions. Species numbers for
bulk and boundary species have to be distinct.
"""
function enable_species!(this::AbstractSystem,ispec::Integer, regions::AbstractVector)
    if ispec>num_species(this)
        throw(DomainError(ispec,"Number of species exceeded"))
    end
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
    function enable_boundary_species(this::AbstractSystem, ispec::Integer, regions::AbstractVector)

Add species to a list of boundary regions. Species numbers for
bulk and boundary species have to be distinct.

"""
function enable_boundary_species!(this::AbstractSystem, ispec::Integer, regions::AbstractVector)
    if ispec>num_species(this)
        throw(DomainError(ispec,"Number of species exceeded"))
    end
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
    isdof(this::AbstractSystem,ispec,inode)
    
Check if degree of freedom is defined.
"""
isdof(this::AbstractSystem,ispec,inode)= this.node_dof[ispec,inode]==ispec ? true : false

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
    num_species(this::AbstractSystem)

Number of species in system
"""
num_species(this::AbstractSystem{Tv}) where Tv = this.physics.num_species



##################################################################
"""
    num_species(this::AbstractSystem)

Retrieve user data record.
"""
data(this::AbstractSystem{Tv}) where Tv = this.physics.data




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










#
# Initialize Dirichlet BC
#
function _inidirichlet!(this::AbstractSystem,U) where Tv
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

#
# Assemble dummy equations for inactive species
#
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

#
# Initialize values in inactive dof for dense system
#
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


#
# Dummy routine for sparse system
#
function _inactspecloop(this::SparseSystem,U,Uold,F)
end

#
# Dummy routine for sparse system
#
function _inactspecinit(this::SparseSystem,U)
end


#
# Accessors for node-dof based loops
#
_firstnodedof(U::SparseSolutionArray{Tv},K) where Tv =U.node_dof.colptr[K]
_lastnodedof(U::SparseSolutionArray{Tv},K) where Tv=U.node_dof.colptr[K+1]-1
_spec(U::SparseSolutionArray{Tv},idof,K) where Tv=U.node_dof.rowval[idof]
_add(U::SparseSolutionArray{Tv},idof,val) where Tv=U.node_dof.nzval[idof]+=val

_firstnodedof(U::Matrix{Tv},K) where Tv = (K-1)*size(U,1)+1
_lastnodedof(U::Matrix{Tv},K) where Tv = K*size(U,1)
_spec(U::Matrix{Tv},idof,K) where Tv =   idof-(K-1)*size(U,1)
_add(U::Matrix{Tv},idof,val) where Tv=U[CartesianIndices(U)[idof]]+=val



##################################################################
"""
    mutable struct BNode

Structure holding local boundary  node information.
Fields:

    index::Int32
    region::Int32
    coord::Array{Tv,1}
    nspec::Int64
"""
mutable struct BNode{Tv}
    index::Int32
    region::Int32
    coord::Array{Tv,1}
    nspec::Int64
    BNode{Tv}(sys::AbstractSystem{Tv}) where Tv  =new(0,0,zeros(Tv,dim_space(sys.grid)),num_species(sys))
end

function _fill!(node::BNode{Tv},grid::Grid{Tv},ibnode,ibface) where Tv
    K=grid.bfacenodes[ibnode,ibface]
    node.region=grid.bfaceregions[ibface]
    node.index=K
    for i=1:length(node.coord)
        node.coord[i]=grid.coord[i,K]
    end
end





##################################################################
"""
    mutable struct Node

Structure holding local node information.
Fields:

    index::Int32
    region::Int32
    coord::Array{Tv,1}
    nspec::Int64
"""
mutable struct Node{Tv}
    index::Int32
    region::Int32
    coord::Array{Tv,1}
    nspec::Int64
    Node{Tv}(sys::AbstractSystem{Tv}) where Tv  =new(0,0,zeros(Tv,dim_space(sys.grid)),num_species(sys))
end

function _fill!(node::Node{Tv},grid::Grid{Tv},inode,icell) where Tv
        K=cellnode(grid,inode,icell)
        node.region=grid.cellregions[icell]
        node.index=K
        for i=1:length(node.coord)
            node.coord[i]=grid.coord[i,K]
        end
end

##################################################################
"""
    mutable struct Edge

Structure holding local edge information.

Fields:

    index::Int32
    nodeK::Int32
    nodeL::Int32
    region::Int32
    coordK::Array{Tv,1}
    coordL::Array{Tv,1}
    nspec::Int64
"""
mutable struct Edge{Tv}
    index::Int32
    nodeK::Int32
    nodeL::Int32
    region::Int32
    coordK::Array{Tv,1}
    coordL::Array{Tv,1}
    nspec::Int64
    Edge{Tv}(sys::AbstractSystem{Tv}) where Tv  =new(0,0,0,0,zeros(Tv,dim_space(sys.grid)),zeros(Tv,dim_space(sys.grid)),num_species(sys))
end


function _fill!(edge::Edge{Tv},grid::Grid{Tv},iedge,icell) where Tv
    K=celledgenode(grid,1,iedge,icell)
    L=celledgenode(grid,2,iedge,icell)
    edge.region=grid.cellregions[icell]
    edge.index=iedge
    edge.nodeK=K
    edge.nodeL=L
    for i=1:length(edge.coordK)
        edge.coordK[i]=grid.coord[i,K]
        edge.coordL[i]=grid.coord[i,L]
    end
end


##################################################################
"""
    num_species(edge::Edge{Tv}) where Tv=edge.nspec

Return number of species for edge
"""
num_species(edge::Edge{Tv}) where Tv=edge.nspec

##################################################################
"""
   function edgelength(edge::Edge)
   
Calculate the length of an edge. 
"""
function edgelength(edge::Edge{Tv}) where Tv
    l::Tv
    l=0.0
    for i=1:length(edge.coordK)
        d=edge.coordK[i]-edge.coordL[i]
        l=l+d*d
    end
    return l
end


@inline viewK(edge::Edge{Tv},u::AbstractArray) where Tv=@views u[1:edge.nspec]
@inline viewL(edge::Edge{Tv},u::AbstractArray) where Tv=@views u[edge.nspec+1:2*edge.nspec]

@inline viewK(nspec::Int64,u::AbstractArray)=@views u[1:nspec]
@inline viewL(nspec::Int64,u::AbstractArray)=@views u[nspec+1:2*nspec]



####################################################################################


