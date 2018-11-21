TwoPointFluxFVM
===============
Solver for coupled nonlinear partial differential equations
based on the two point flux finite volume method on admissible grids.


This Julia package merges the ideas behind [pdelib](http://www.wias-berlin.de/software/pdelib/?lang=0)/fvsys with the multiple dispatch paradigm of Julia. It allows to avoid the implementation of function derivatives. It instead uses the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and [DiffResults](https://github.com/JuliaDiff/DiffResults.jl)  to evaluate user fucntions along with their jacobians.

So far this is merely a design study to learn and evaluate Julia.  
It is however aimed to be feasible at least for small projects.

It requires Julia 1.0.

Documentation created with [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/index.html)
resides [here](https://www.wias-berlin.de/people/fuhrmann/TwoPointFluxFVM)

## Main structs (classes)

### [`TwoPointFluxFVM.Physics`](@ref)
This is an abstract type  from which a user
data type can be derived which describes the physical
data of a given problem. These are user defined parameters
and nonlinear functions describing storage, reactions and
fluxes between control volumes.

A typical usage pattern is as follows:

```julia
mutable struct Physics <:TwoPointFluxFVM.Physics
    TwoPointFluxFVM.@AddPhysicsBaseClassFields
    eps::Float64 
    Physics()=Physics(new())
end

function reaction!(this::Physics,f::AbstractArray,u::AbstractArray)
    f[1]=u[1]^2
end

function flux!(this::Physics,f::AbstractArray,uk::AbstractArray,ul::AbstractArray)
    f[1]=this.eps*(uk[1]^2-ul[1]^2)
end 

function source!(this::Physics,f::AbstractArray,x::AbstractArray)
    f[1]=1.0e-4*x[1]
end 

function Physics(this::Physics)
    TwoPointFluxFVM.PhysicsBase(this,1)
    this.eps=1
    this.flux=flux!
    this.reaction=reaction!
    this.source=source!
    return this
end
```

### [`TwoPointFluxFVM.Graph`](@ref)

This is a weighted graph which represents edges and nodes
of a finite volume scheme. There are currently two
constructors:

[`TwoPointFluxFVM.Graph(X::Array{Float64,1})`](@ref) for one-dimensional
domains and [`TwoPointFluxFVM.Graph(X::Array{Float64,1},Y::Array{Float64,1})`](@ref)
for two-dimensional domains.

### [`TwoPointFluxFVM.System`](@ref)

From instances of  [`TwoPointFluxFVM.Graph`](@ref) and [`TwoPointFluxFVM.Physics`](@ref),
a [`TwoPointFluxFVM.System`](@ref) which contains all the necessary
data for the solution of the nonlinear system described by them.




