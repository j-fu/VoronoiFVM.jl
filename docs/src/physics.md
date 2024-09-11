# Physics & special functions

## Physics
```@docs
VoronoiFVM.AbstractPhysics
VoronoiFVM.Physics
VoronoiFVM.Physics(;kwargs...)

Base.show(io::IO,physics::VoronoiFVM.AbstractPhysics)
```
## Handling boundary conditions
Boundary conditions are handled in the  `bcondition` callback passed to the system constructor.
For being called in this callback, the following  functions are available

```@docs
boundary_dirichlet!(y,u,bnode::AbstractGeometryItem,ispec,ireg,val)
boundary_dirichlet!(y,u,bnode::AbstractGeometryItem;kwargs...)
boundary_neumann!(y,u,bnode::AbstractGeometryItem,ispec,ireg,val)
boundary_neumann!(y,u,bnode::AbstractGeometryItem;kwargs...)
boundary_robin!(y,u,bnode::AbstractGeometryItem,ispec,ireg,fac,val)
boundary_robin!(y,u,bnode::AbstractGeometryItem;kwargs...)
ramp
```

### Outflow boundary conditions
These are characterized by the `boutflow` physics callback and 
and the `outflowboundaries` keyword argument in the system
resp. physics constructor. See also the 
[corresponding notebook](https://j-fu.github.io/VoronoiFVM.jl/dev/nbhtml/outflow/)

```@docs
hasoutflownode
isoutflownode
outflownode
```

## Edge and node data
```@docs
VoronoiFVM.Node
VoronoiFVM.BNode
VoronoiFVM.Edge
VoronoiFVM.BEdge
```

## Special functions
```@docs
fbernoulli
fbernoulli_pm
inplace_linsolve!(A,b)
inplace_linsolve!(A,b,ipiv)
```
