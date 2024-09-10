"""
```
const DenseSolutionArray=Matrix
```

Dense storage of solution
"""
mutable struct DenseSolutionArray{T, N}  <: AbstractSolutionArray{T,N}
    u::Array{T,N}
    history::Union{NewtonSolverHistory, Nothing}
end

DenseSolutionArray(u::Matrix{T}) where {T} =DenseSolutionArray{T,2}(u,nothing)

Base.getindex(a::DenseSolutionArray, i::Int, j::Int)= getindex(a.u,i,j )
Base.setindex!(a::DenseSolutionArray,v, i::Int, j::Int) = setindex!(a.u,v,i,j)
Base.size(a::DenseSolutionArray)=size(a.u)

"""
 $(SIGNATURES)
    
Get degree of freedom number
"""
dof(a::DenseSolutionArray, ispec, K) = (K - 1) * size(a, 1) + ispec

"""
$(SIGNATURES)

Return indices for dense solution array.
"""
unknown_indices(a::DenseSolutionArray) = LinearIndices(a)



Base.vec(a::DenseSolutionArray) = vec(a.u)

Base.copy(a::DenseSolutionArray) = DenseSolutionArray(copy(a.u))

Base.similar(a::DenseSolutionArray) = DenseSolutionArray(similar(a.u))

"""
$(SIGNATURES)

Array of values in solution array.
"""
values(a::DenseSolutionArray) = vec(a)

"""
$(TYPEDSIGNATURES)

Add residual value into global degree of freedom
"""
_add(U::DenseSolutionArray, idof, val) = U[CartesianIndices(U)[idof]] += val

_set(U::DenseSolutionArray, idof, val) = U[CartesianIndices(U)[idof]] = val

##################################################################
