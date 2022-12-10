"""
```
const DenseSolutionArray=Matrix
```

Dense storage of solution
"""
const DenseSolutionArray=Matrix


"""
$(SIGNATURES)

Get degree of freedom number
"""
dof(a::DenseSolutionArray{Tv}, ispec::Integer, K::Integer) where Tv = (K-1)*size(a,1)+ispec

unknown_indices(a::DenseSolutionArray{Tv}) where Tv=LinearIndices(a)


"""
$(SIGNATURES)

Array of values in solution array.
"""
values(a::DenseSolutionArray{Tv}) where Tv = vec(a)


#
# Accessors for node-dof based loops
#
_firstnodedof(U::DenseSolutionArray{Tv},K::Integer) where Tv = (K-1)*size(U,1)+1
_lastnodedof(U::DenseSolutionArray{Tv},K::Integer) where Tv = K*size(U,1)
_species_of_dof(U::DenseSolutionArray{Tv},idof,K) where Tv =   idof-(K-1)*size(U,1)
_add(U::DenseSolutionArray{Tv},idof,val) where Tv=U[CartesianIndices(U)[idof]]+=val


##################################################################

