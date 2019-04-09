VoronoiFVM
===============
Solver for coupled nonlinear partial differential equations
based on the two point flux finite volume method on admissible grids.


This Julia package merges the ideas behind [pdelib](http://www.wias-berlin.de/software/pdelib/?lang=0)/fvsys with the multiple dispatch paradigm of Julia. It allows to avoid the implementation of function derivatives. It instead uses the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and [DiffResults](https://github.com/JuliaDiff/DiffResults.jl)  to evaluate user fucntions along with their jacobians.

So far this is merely a design study to learn and evaluate Julia.  
It is however aimed to be feasible at least for small projects.

It requires Julia 1.0.

Documentation created with [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/index.html)
resides [here](https://www.wias-berlin.de/people/fuhrmann/VoronoiFVM)

## Typical usage

```julia
"""
Structure containing  userdata information
"""
mutable struct Physics
    reaction::Function
    flux::Function
    storage::Function
    source::Function
    eps::Float64 
    Physics()=new()
end

"""
Reaction term
"""
physics.reaction=function(physics,node,f,u)
    f[1]=u[1]^2
end

"""
Storage term
"""
physics.storage=function(physics,node,f,u)
    f[1]=u[1]
end

"""
Flux term
"""
physics.flux=function(physics,edge,f,uk,ul)
    f[1]=this.eps*(uk[1]^2-ul[1]^2)
end 


"""
Source term
"""
physics.source=function(physics,node,f)
    f[1]=1.0e-4*node.coord[1]
end 


```

### [`VoronoiFVM.Grid`](@ref)

This is a simplex grid structure.

### [`VoronoiFVM.Node`](@ref)

This represents a node  in [`VoronoiFVM.Grid`](@ref).

### [`VoronoiFVM.Edge`](@ref)

This represents an edge between
two neigboring control volumes created from [`VoronoiFVM.Grid`](@ref).

Currently, constructors are
[`VoronoiFVM.Grid(X::Array{Real,1})`](@ref) for one-dimensional
domains and [`VoronoiFVM.Grid(X::Array{Float64,1},Y::Array{Float64,1})`](@ref)
for two-dimensional domains.

### [`VoronoiFVM.System`](@ref)

From instances of  [`VoronoiFVM.Graph`](@ref) and [`VoronoiFVM.Physics`](@ref),
a [`VoronoiFVM.System`](@ref) which contains all the necessary
data for the solution of the nonlinear system described by them.




