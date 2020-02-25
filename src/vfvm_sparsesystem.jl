##################################################################
"""
$(TYPEDEF)

Struct holding solution information for SparseSystem. Solution
is stored in a sparse matrix structure.

This class plays well with the abstract array interface.

$(TYPEDFIELDS)
"""
struct SparseSolutionArray{Tv,Ti} <: AbstractMatrix{Tv}

    """
    Sparse matrix holding actual data.
    """
    node_dof::SparseMatrixCSC{Tv,Ti}
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
of the sparse matrix structures for the bookkeeping of the unknowns
creates overhead.

$(TYPEDFIELDS)
"""
mutable struct SparseSystem{Tv,Ti} <: AbstractSystem{Tv,Ti}

    """
    Grid
    """
    grid::Grid{Tv, Ti}

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
    region_species::SparseMatrixCSC{Int8,Ti}

    """
    Sparse matrix containing species numbers for boundary regions
    """
    bregion_species::SparseMatrixCSC{Int8,Ti}


    """
    Sparse matrix containing degree of freedom numbers for each node
    """
    node_dof::SparseMatrixCSC{Int8,Ti}

    """
    Jacobi matrix for nonlinear problem
    """
    matrix::ExtendableSparseMatrix{Tv,Ti}

    """
    Flag which says if the number of unknowns per node is constant
    """
    species_homogeneous::Bool

    """
    Solution vector holding Newton update
    """
    update::SparseSolutionArray{Tv,Ti}

    """
    Solution vector holding Newton residual
    """
    residual::SparseSolutionArray{Tv,Ti}

    """
    API version
    """
    oldapi::Bool

    SparseSystem{Tv,Ti}() where {Tv,Ti} = new()
end

##################################################################
"""
$(TYPEDSIGNATURES)

Constructor for SparseSystem. `physics` provides some user data, `maxspec`
is the maximum number of species.
"""
function  SparseSystem(grid::Grid,physics::Physics; oldapi=true)
    Tv=Base.eltype(grid)
    Ti=eltype(grid.cellnodes)
    this=SparseSystem{Tv,Ti}()
    maxspec=physics.num_species
    this.grid=grid
    this.physics=physics
    this.region_species=spzeros(Int8,Int16,maxspec,num_cellregions(grid))
    this.bregion_species=spzeros(Int8,Int16,maxspec,num_bfaceregions(grid))
    this.node_dof=spzeros(Int8,Int32,maxspec,num_nodes(grid))
    this.boundary_values=zeros(Tv,maxspec,num_bfaceregions(grid))
    this.boundary_factors=zeros(Tv,maxspec,num_bfaceregions(grid))
    this.species_homogeneous=false
    this.oldapi=oldapi
    return this
end


##################################################################
"""
$(TYPEDSIGNATURES)

Number of degrees of freedom for system.
"""
num_dof(this::SparseSystem)= nnz(this.node_dof)



##################################################################
"""
$(SIGNATURES)

Create a solution vector for system. 
The entries of the returned vector are undefined.
"""
function unknowns(sys::SparseSystem{Tv,Ti}) where {Tv,Ti}
    return SparseSolutionArray{Tv,Ti}(SparseMatrixCSC(sys.node_dof.m,
                                                   sys.node_dof.n,
                                                   sys.node_dof.colptr,
                                                   sys.node_dof.rowval,
                                                   Array{Tv}(undef,num_dof(sys))
                                                   )
                                   )
end

##################################################################
"""
$(SIGNATURES)

Create a solution vector for system with given type
The entries of the returned vector are undefined.
"""
function unknowns(Tu::Type, sys::SparseSystem{Tv, Ti}) where {Tv,Ti}
    return SparseSolutionArray{Tu,Ti}(SparseMatrixCSC(sys.node_dof.m,
                                                      sys.node_dof.n,
                                                      sys.node_dof.colptr,
                                                      sys.node_dof.rowval,
                                                      Array{Tu}(undef,num_dof(sys))
                                                      )
                                   )
end



##################################################################
"""
$(SIGNATURES)

Reshape vector to fit as solution to system.
"""
function Base.reshape(v::AbstractVector{Tu},sys::SparseSystem{Tv,Ti}) where {Tu,Tv,Ti}
    @assert  length(v)==num_dof(sys)
    SparseSolutionArray{Tu,Ti}(SparseMatrixCSC(sys.node_dof.m,
                                            sys.node_dof.n,
                                            sys.node_dof.colptr,
                                            sys.node_dof.rowval,
                                            Vector{Tu}(v)
                                            )
                            )
end

Base.reshape(v::SparseSolutionArray,sys::SparseSystem)=v



##################################################################
"""
$(TYPEDSIGNATURES)
    
Return size of solution array.
"""
Base.size(a::SparseSolutionArray)=size(a.node_dof)


##################################################################
"""
$(TYPEDSIGNATURES)

Array of values in solution array.
"""
values(a::SparseSolutionArray)=a.node_dof.nzval




##################################################################
"""
$(SIGNATURES)

Create a copy of solution array
"""
Base.copy(this::SparseSolutionArray{Tv,Ti}) where {Tv,Ti} = SparseSolutionArray{Tv,Ti}(SparseMatrixCSC(this.node_dof.m,
                                                                                                       this.node_dof.n,
                                                                                                       this.node_dof.colptr,
                                                                                                       this.node_dof.rowval,
                                                                                                       Base.copy(this.node_dof.nzval)
                                                                                                       )
                                                                            )
##################################################################
"""
$(SIGNATURES)

Get number of degree of freedom. Return 0 if species is not defined in node.
"""
@inline function dof(a::SparseSolutionArray{Tv,Ti},i::Integer, j::Integer) where {Tv,Ti}
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
_firstnodedof(U::SparseSolutionArray,K) =U.node_dof.colptr[K]
_lastnodedof(U::SparseSolutionArray,K) =U.node_dof.colptr[K+1]-1
_spec(U::SparseSolutionArray,idof,K) =U.node_dof.rowval[idof]
_add(U::SparseSolutionArray,idof,val)=U.node_dof.nzval[idof]+=val


