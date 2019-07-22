##########################################################
"""
$(TYPEDEF)
    
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
$(TYPEDEF)

Struct holding solution information for SparseSystem. Solution
is stored in a sparse matrix structure.

This class plays well with the abstract array interface.

$(TYPEDFIELDS)
"""
struct SparseSolutionArray{Tv} <: AbstractMatrix{Tv}

    """
    Sparse matrix holding actual data.
    """
    node_dof::SparseMatrixCSC{Tv,Int32}
end


##################################################################
"""
$(TYPEDEF)

Structure holding data for finite volume system solution.
Information on species distribution is kept in sparse
matrices, and the solution array is of type SparseSolutionArray,
i.e. effectively it is a sparse matrix.

Unlike in the DenseSystem, the system matrix handles exactly those
degrees of freedom which correspond to unknowns. However, handling
of the sparse matrix structures for the bookeeping of the unknowns
creates overhead.

$(TYPEDFIELDS)
"""
mutable struct SparseSystem{Tv} <: AbstractSystem{Tv}

    """
    Grid
    """
    grid::Grid

    """
    Physics data
    """
    physics::Physics

    """
    Array of boundary values 
    """
    boundary_values::Array{Tv,2} 

    """
    Array of boundary factors 
    """
    boundary_factors::Array{Tv,2}

    """
    Sparse matrix containing species numbers for inner regions
    """
    region_species::SparseMatrixCSC{Int8,Int16}

    """
    Sparse matrix containing species numbers for boundary regions
    """
    bregion_species::SparseMatrixCSC{Int8,Int16}


    """
    Sparse matrix containing degree of freedom numbers for each node
    """
    node_dof::SparseMatrixCSC{Int8,Int32}

    """
    Jacobi matrix for nonlinear problem
    """
    matrix::ExtendableSparseMatrix{Tv,Int64}

    """
    Flag which says if the number of unknowns per node is constant
    """
    species_homogeneous::Bool

    """
    Solution vector holding Newton update
    """
    update::SparseSolutionArray{Tv}

    """
    Solution vector holding Newton residual
    """
    residual::SparseSolutionArray{Tv}

    SparseSystem{Tv}() where Tv = new()
end

##################################################################
"""
$(TYPEDSIGNATURES)

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
$(TYPEDEF)

Structure holding data for finite volume system solution.
Information on species distribution is kept in dense
matrices, and the solution array is of type Array{2}.


Unlike in the SparseSystem, the system matrix handles exactly those
degrees of freedom which correspond to unknowns, and dummy 
degrees of freedom where unknowns are not defined. Handling
of the sparse matrix structures for the bookeeping of the unknowns
has less overhead, but additional dummy equations are added
to the system matrix.

$(TYPEDFIELDS)
"""
mutable struct DenseSystem{Tv} <: AbstractSystem{Tv}
    """
    Grid
    """
    grid::Grid

    """
    Physics data
    """
    physics::Physics

    """
    Array of boundary values 
    """
    boundary_values::Array{Tv,2} 

    """
    Array of boundary factors
    """
    boundary_factors::Array{Tv,2}

    """
    Full matrix containing species numbers for inner regions
    """
    region_species::Array{Int8,2}

    """
    Full matrix containing species numbers for boundary regions
    """
    bregion_species::Array{Int8,2}

    """
    Full matrix containing degree of freedom numbers for each node
    """
    node_dof::Array{Int8,2}

    """
    Jacobi matrix for nonlinear problem
    """
    matrix::ExtendableSparseMatrix{Tv,Int64}

    """
    Flag which says if the number of unknowns per node is constant
    """
    species_homogeneous::Bool

    """
    Solution vector holding Newton update
    """
    update::Matrix{Tv}

    """
    Solution vector holding Newton residual
    """
    residual::Matrix{Tv}

    DenseSystem{Tv}() where Tv = new()
end
##################################################################
"""
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

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


# Create matrix in system and figure out if species
# distribution is homgeneous
function _complete!(this::AbstractSystem{Tv};create_newtonvectors=false) where Tv
    if isdefined(this,:matrix)
        return
    end
    this.matrix=ExtendableSparseMatrix{Tv,Int64}(num_dof(this), num_dof(this))
    this.species_homogeneous=true
    species_added=false
    for inode=1:size(this.node_dof,2)
        for ispec=1:size(this.node_dof,1)
            if this.node_dof[ispec,inode]==ispec
                species_added=true
            else
                this.species_homogeneous=false
            end
        end
    end
    
    if (!species_added)
        error("No species enabled.\n Call enable_species(system,species_number, list_of_regions) at least once.")
    end
    if create_newtonvectors
        this.residual=unknowns(this)
        this.update=unknowns(this)
    end
end

##################################################################
"""
$(TYPEDSIGNATURES)
    
Check if degree of freedom is defined.
"""
isdof(this::AbstractSystem,ispec,inode)= this.node_dof[ispec,inode]==ispec ? true : false

##################################################################
"""
$(TYPEDSIGNATURES)

Number of degrees of freedom for system.
"""
num_dof(this::SparseSystem)= nnz(this.node_dof)

##################################################################
"""
$(TYPEDSIGNATURES)

Number of degrees of freedom for system.
"""
num_dof(this::DenseSystem)= length(this.node_dof)

##################################################################
"""
$(TYPEDSIGNATURES)

Number of species in system
"""
num_species(this::AbstractSystem{Tv}) where Tv = this.physics.num_species



##################################################################
"""
$(TYPEDSIGNATURES)

Retrieve user data record.
"""
data(this::AbstractSystem{Tv}) where Tv = this.physics.data




##################################################################
"""
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

Create a solution vector for system.
"""
unknowns(sys::DenseSystem{Tv}) where Tv=Array{Tv}(undef,size(sys.node_dof,1), size(sys.node_dof,2))


##################################################################
"""
$(TYPEDSIGNATURES)
    
Return size of solution array.
"""
Base.size(a::SparseSolutionArray)=size(a.node_dof)

##################################################################
"""
$(TYPEDSIGNATURES)
                        
Number of nodes (size of second dimension) of solution array.
"""
num_nodes(a)=size(a,2)

##################################################################
"""
$(TYPEDSIGNATURES)

Number of species (size of first dimension) of solution array.
"""
num_species(a)=size(a,1)

##################################################################
"""
$(TYPEDSIGNATURES)

Array of values in solution array.
"""
values(a::SparseSolutionArray)=a.node_dof.nzval


##################################################################
"""
$(TYPEDSIGNATURES)

Array of values in solution array.
"""
values(a::Array)= vec(a)



##################################################################
"""
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

Get number of degree of freedom.
"""
dof(a::Array{Tv,2}, ispec::Integer, K::Integer) where Tv = (K-1)*size(a,1)+ispec


##################################################################
"""
$(TYPEDSIGNATURES)

Set value for degree of freedom.
"""
function setdof!(a::SparseSolutionArray,v,i::Integer)
    a.node_dof.nzval[i] = v
end

##################################################################
"""
$(TYPEDSIGNATURES)

Return  value for degree of freedom.
"""
getdof(a::SparseSolutionArray,i::Integer) =a.node_dof.nzval[i] 




##################################################################
"""
$(TYPEDSIGNATURES)

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
$(TYPEDSIGNATURES)

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
$(TYPEDEF)

Struct holding information for solution array view on subgrid

$(TYPEDFIELDS)
"""
struct SubgridArrayView{Tv,Ta} <: AbstractMatrix{Tv}

    """
    Original array
    """
    sysarray::Ta

    """
    Subgrid for view
    """
    subgrid::SubGrid
end

##################################################################
"""
$(TYPEDSIGNATURES)

Create a view of the solution array on a subgrid.
"""
Base.view(a::AbstractMatrix{Tv},sg::SubGrid) where Tv = SubgridArrayView{Tv,typeof(a)}(a,sg)


##############################################################################
"""
$(TYPEDSIGNATURES)

Accessor method for subgrid array view.
"""
Base.getindex(aview::SubgridArrayView,ispec::Integer,inode::Integer) = aview.sysarray[ispec,aview.subgrid.node_in_parent[inode]]

##############################################################################
"""
$(TYPEDSIGNATURES)

Accessor method for subgrid array view.
"""
@inline function Base.setindex!(aview::SubgridArrayView,v,ispec::Integer,inode::Integer)
    aview.sysarray[ispec,aview.subgrid.node_in_parent[inode]]=v
    return aview
end

##################################################################
"""
$(TYPEDSIGNATURES)
    
Return size of solution array view.
"""
Base.size(a::SubgridArrayView)=(size(a.sysarray,1),size(a.subgrid.node_in_parent,1))










#
# Initialize Dirichlet BC
#
function _initialize_dirichlet!(U::AbstractMatrix{Tv},this::AbstractSystem{Tv}) where {Tv}
    for ibface=1:num_bfaces(this.grid)
        ibreg=this.grid.bfaceregions[ibface]
        for ispec=1:num_species(this)
            if this.boundary_factors[ispec,ibreg]≈ Dirichlet
                for inode=1:dim_grid(this.grid)
                    U[ispec,this.grid.bfacenodes[inode,ibface]]=this.boundary_values[ispec,ibreg]
                end
            end
        end
    end
end



#
# Assemble dummy equations for inactive species
#
function _eval_and_assemble_inactive_species(this::DenseSystem,U,Uold,F)
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
function _initialize_inactive_dof!(U::Matrix{Tv},this::DenseSystem{Tv}) where {Tv}
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

function _initialize!(U::AbstractMatrix{Tv},this::AbstractSystem{Tv}) where {Tv}
    _initialize_dirichlet!(U,this)
    _initialize_inactive_dof!(U,this)
end


#
# Dummy routine for sparse system
#
function _eval_and_assemble_inactive_species(this::SparseSystem,U,Uold,F)
end

#
# Dummy routine for sparse system
#
function     _initialize_inactive_dof!(U::AbstractMatrix{Tv},this::SparseSystem{Tv}) where {Tv}
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
$(TYPEDEF)

Structure holding local boundary  node information.

$(TYPEDFIELDS)
"""
mutable struct BNode{Tv}

    """
    Index in grid
    """
    index::Int32

    """
    Boundary region number
    """
    region::Int32

    """
    1D Array of node coordinates
    """
    coord::Array{Tv,1}

    """
    Number of species defined in node
    """
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
$(TYPEDEF)

Structure holding local node information.

$(TYPEDFIELDS)
"""
mutable struct Node{Tv}

    """
    Index in grid

    """
    index::Int32
    """
    Inner region number
    """
    region::Int32

    """
    1D Array of node coordinates
    """
    coord::Array{Tv,1}

    """
    Number of species defined in node
    """
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
$(TYPEDEF)

Structure holding local edge information.

$(TYPEDFIELDS)
"""
mutable struct Edge{Tv}

    """
    Index in grid
    """
    index::Int32

    """
    Index of first node
    """
    nodeK::Int32

    """
    Index of second node
    """
    nodeL::Int32

    """
    Inner region number corresponding to edge
    """
    region::Int32

    """
    1D Array of first node coordinates
    """
    coordK::Array{Tv,1}

    """
    1D Array of second node coordinates
    """
    coordL::Array{Tv,1}

    """
    Number of species defined in edge
    """
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
$(TYPEDSIGNATURES)

Return number of species for edge
"""
num_species(edge::Edge{Tv}) where Tv=edge.nspec

##################################################################
"""
$(TYPEDSIGNATURES)
   
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

"""
$(TYPEDSIGNATURES)

Solution view on first edge node
"""
@inline viewK(edge::Edge{Tv},u::AbstractArray) where Tv=@views u[1:edge.nspec]


"""
$(TYPEDSIGNATURES)

Solution view on second edge node
"""
@inline viewL(edge::Edge{Tv},u::AbstractArray) where Tv=@views u[edge.nspec+1:2*edge.nspec]

"""
$(TYPEDSIGNATURES)

Solution view on first edge node
"""
@inline viewK(nspec::Int64,u::AbstractArray)=@views u[1:nspec]


"""
$(TYPEDSIGNATURES)

Solution view on second edge node
"""
@inline viewL(nspec::Int64,u::AbstractArray)=@views u[nspec+1:2*nspec]



function Base.show(io::IO,sys::AbstractSystem) where Tc
    str=@sprintf("%s(num_species=%d)",typeof(sys),sys.physics.num_species)
    println(io,str)
end


####################################################################################


