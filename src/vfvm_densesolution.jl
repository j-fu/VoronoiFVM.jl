"""
```
const DenseSolutionArray=Matrix
```

Dense storage of solution
"""
const DenseSolutionArray = Matrix

"""
$(SIGNATURES)

Get degree of freedom number
"""
dof(a::DenseSolutionArray{Tv}, ispec::Integer, K::Integer) where {Tv} = (K - 1) * size(a, 1) + ispec

unknown_indices(a::DenseSolutionArray{Tv}) where {Tv} = LinearIndices(a)

"""
$(SIGNATURES)

Array of values in solution array.
"""
values(a::DenseSolutionArray{Tv}) where {Tv} = vec(a)

"""
$(TYPEDSIGNATURES)

Add residual value into global degree of freedom
"""
_add(U::DenseSolutionArray{Tv}, idof, val) where {Tv} = U[CartesianIndices(U)[idof]] += val

_set(U::DenseSolutionArray{Tv}, idof, val) where {Tv} = U[CartesianIndices(U)[idof]] = val

##################################################################
