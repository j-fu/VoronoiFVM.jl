
#################################################################
"""
$(TYPEDEF)

Struct holding solution information for SparseSystem. Solution
is stored in a sparse matrix structure.

This class plays well with the abstract array interface.

Fields:

$(TYPEDFIELDS)
"""
struct SparseSolutionArray{Tv, Ti} <: AbstractMatrix{Tv}
    """
    Sparse matrix holding actual data.
    """
    node_dof::SparseMatrixCSC{Tv, Ti}
end

##################################################################
"""
$(SIGNATURES)
    
Return size of sparse solution array.
"""
Base.size(a::SparseSolutionArray) = size(a.node_dof)


##################################################################
"""
$(SIGNATURES)

Array of values in sparse solution array.
"""
values(a::SparseSolutionArray) = a.node_dof.nzval

##################################################################
"""
$(SIGNATURES)
    
Create a copy of sparse solution array
"""
function Base.copy(this::SparseSolutionArray{Tv, Ti}) where {Tv, Ti}
    SparseSolutionArray{Tv, Ti}(SparseMatrixCSC(this.node_dof.m,
                                                this.node_dof.n,
                                                this.node_dof.colptr,
                                                this.node_dof.rowval,
                                                Base.copy(this.node_dof.nzval)))
end

"""
$(SIGNATURES)
    
Create a similar uninitialized sparse solution array
"""
function Base.similar(this::SparseSolutionArray{Tv, Ti}) where {Tv, Ti}
    SparseSolutionArray{Tv, Ti}(SparseMatrixCSC(this.node_dof.m,
                                                this.node_dof.n,
                                                this.node_dof.colptr,
                                                this.node_dof.rowval,
                                                Base.similar(this.node_dof.nzval)))
end
##################################################################
"""
$(SIGNATURES)

Get number of degree of freedom. Return 0 if species is not defined in node.
"""
function dof(a::SparseSolutionArray{Tv, Ti}, i::Integer, j::Integer) where {Tv, Ti}
    A = a.node_dof
    coljfirstk = Int(A.colptr[j])
    coljlastk = Int(A.colptr[j + 1] - 1)
    searchk = searchsortedfirst(A.rowval, i, coljfirstk, coljlastk, Base.Order.Forward)
    if searchk <= coljlastk && A.rowval[searchk] == i
        return searchk
    end
    return 0
end

struct SparseSolutionIndices
    a::SparseSolutionArray
end

unknown_indices(a::SparseSolutionArray) = SparseSolutionIndices(a)

Base.getindex(idx::SparseSolutionIndices, i, j) = dof(idx.a, i, j)

##################################################################
"""
$(SIGNATURES)

Set value for degree of freedom.
"""
function setdof!(a::SparseSolutionArray, v, i::Integer)
    a.node_dof.nzval[i] = v
end

##################################################################
"""
$(SIGNATURES)

Return  value for degree of freedom.
"""
getdof(a::SparseSolutionArray, i::Integer) = a.node_dof.nzval[i]

Base.:-(a::SparseSolutionArray, b::SparseSolutionArray) = SparseSolutionArray(a.node_dof - b.node_dof)
Base.:+(a::SparseSolutionArray, b::SparseSolutionArray) = SparseSolutionArray(a.node_dof + b.node_dof)

##################################################################
"""
$(SIGNATURES)

Accessor for sparse solution array.
"""
function Base.setindex!(a::SparseSolutionArray, v, ispec::Integer, inode::Integer)
    searchk = dof(a, ispec, inode)
    if searchk > 0
        setdof!(a, v, searchk)
        return a
    end
    # TODO: what is the right reacton here ?
    # Ignoring seems to be better, so we can broacast etc.
    # throw(DomainError("undefined degree of freedom"))
end

##################################################################
"""
$(SIGNATURES)

Accessor for sparse solution array.
"""
function Base.getindex(a::SparseSolutionArray, ispec::Integer, inode::Integer)
    searchk = dof(a, ispec, inode)
    if searchk > 0
        return getdof(a, searchk)
    end
    #
    # TODO: what is the right reacton here ?
    # Actually, NaN plays well with pyplot...
    return NaN
end

"""
$(TYPEDSIGNATURES)

Add residual value into global degree of freedom
"""
_add(U::SparseSolutionArray, idof, val) = U.node_dof.nzval[idof] += val

_set(U::SparseSolutionArray, idof, val) = U.node_dof.nzval[idof] = val
