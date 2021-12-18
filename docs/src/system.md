# System

The computational grid required is assumed to correspond to a domain
``\Omega=\cup_{r=1}^{n_\Omega} \Omega_r`` 

Grids for VoronoiFVM are managed by the packages [ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl)
and [SimplexGridFactory.jl](https://github.com/j-fu/SimplexGridFactory.jl)


with boundary  ``\partial\Omega=\Gamma=\cup_{b=1}^{n_\Gamma} \Gamma_b``.

The subdomains ``\Omega_r`` are called "regions" and the boundary subdomains ``\Gamma_b`` are called "boundary regions".

On this complex of domains "lives"  a number of species which are either attached to a number of regions or to a number of boundary regions.

All these data, the matrix for the linear system and other things are hold together by a struct `VoronoiFVM.System`. 
This type is not exported to avoid name clashes.


## System constructors

```@docs
VoronoiFVM.System(grid::ExtendableGrid; kwargs...)
VoronoiFVM.System(X::AbstractVector; kwargs...)
VoronoiFVM.System(X::AbstractVector,Y::AbstractVector; kwargs...)
VoronoiFVM.System(X::AbstractVector,Y::AbstractVector,Z::AbstractVector; kwargs...)
update_grid!
```

## Adding species by species numbers
```@docs
enable_species!(system::VoronoiFVM.AbstractSystem,ispec::Integer, regions::AbstractVector)
enable_species!(system::VoronoiFVM.AbstractSystem; kwargs...)
enable_boundary_species!
```


## Handling boundary conditions
Boundary conditions are handeled in the  `bcondition` callback passed to the system constructor.
For being called in this callback, the following  functions are availabel

```@docs
boundary_dirichlet!(y,u,bnode,ispec,ireg,val)
boundary_dirichlet!(y,u,bnode;kwargs...)
boundary_neumann!(y,u,bnode,ispec,ireg,val)
boundary_neumann!(y,u,bnode;kwargs...)
boundary_robin!(y,u,bnode,ispec,ireg,fac,val)
boundary_robin!(y,u,bnode;kwargs...)
ramp
```

## Various tools

```@docs
physics!
```

```@docs
check_allocs!
```


```@docs
num_dof
unknowns(system::VoronoiFVM.AbstractSystem; kwargs...)
VoronoiFVM.isunknownsof
Base.reshape
LinearAlgebra.norm(system::VoronoiFVM.AbstractSystem, u,p)
```

## Types

```@docs
VoronoiFVM.AbstractSystem
VoronoiFVM.System{Tv,Ti, Tm, TSpecMat<:AbstractMatrix, TSolArray<:AbstractMatrix}
VoronoiFVM.DenseSystem
VoronoiFVM.SparseSystem
```


## Legacy API
```@docs
boundary_dirichlet!(system::VoronoiFVM.AbstractSystem, ispec, ibc, v)
boundary_dirichlet!(system::VoronoiFVM.AbstractSystem; kwargs...)
boundary_neumann!(system::VoronoiFVM.AbstractSystem, ispec, ibc, v)
boundary_neumann!(system::VoronoiFVM.AbstractSystem; kwargs...)
boundary_robin!(system::VoronoiFVM.AbstractSystem, ispec, ibc,alpha, v)
boundary_robin!(system::VoronoiFVM.AbstractSystem; kwargs...)
```

