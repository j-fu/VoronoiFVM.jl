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
physics!
```

## Adding species by species numbers
```@docs
enable_species!(system::VoronoiFVM.AbstractSystem,ispec::Integer, regions::AbstractVector)
enable_species!(system::VoronoiFVM.AbstractSystem; kwargs...)
enable_boundary_species!
```


## Handling boundary conditions
Boundary conditions are handled in the  `bcondition` callback passed to the system constructor.
For being called in this callback, the following  functions are available

```@docs
boundary_dirichlet!(y,u,bnode,ispec,ireg,val)
boundary_dirichlet!(y,u,bnode;kwargs...)
boundary_neumann!(y,u,bnode,ispec,ireg,val)
boundary_neumann!(y,u,bnode;kwargs...)
boundary_robin!(y,u,bnode,ispec,ireg,fac,val)
boundary_robin!(y,u,bnode;kwargs...)
ramp
```

## Allocation warnings

The code checks for allocations in the assembly loop. 
Care has been taken to ensure that allocations in the assembly loop don't emerge
from VoronoiFVM.jl code.

If allocations occur in the assembly  loop, they happen in the physics
callbacks.  The corresponding warnings can bee switched off by passing
a  verbosity strings  without  'a'  to the  solver.   If  no data  are
allocated in the physics callbacks, these allocations are probably due to 
type instabilities in physics callbacks, see the the discussion
[here](../runexamples/#Performance-with-closures).  Type instabilities
can be debugged via the `@time`  macro applied to expressions in a
physics callback.

The following  cases provide some ideas  where to look for  reasons of
the problem and possible remedies:

Case 1: a parameter changes its value, and Julia is not sure about the type.
```julia
eps=1.0

flux(f,_u,edge)
    u=unknowns(edge,_u)
    f[1]=eps*(u[1,1]-[1,2])
end
... solve etc ...
eps=2.0
```
Remedy: use a type annotation `eps::Float64=...` to signalize your intent to Julia.
This behaviour is explained in the [Julia documentation](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured).



Case 2: variables in the closure have the same name as a variable
introduced in a callback.
```julia
flux(f,_u,edge)
    u=unknowns(edge,_u)
    f[1]=(u[1,1]-[1,2])
end

... create etc ...

u=solve(...)
```
Remedy: rename e.g. `u=solve()` to `sol=solve()`



## Various tools

```@docs
num_dof
unknowns(system::VoronoiFVM.AbstractSystem; kwargs...)
Base.map
Base.map!
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

