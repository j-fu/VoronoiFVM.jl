# Physics & special functions

## Physics
```@docs
VoronoiFVM.AbstractPhysics
VoronoiFVM.Physics
VoronoiFVM.Physics(;kwargs...)
Base.show(io::IO,physics::VoronoiFVM.AbstractPhysics)
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
