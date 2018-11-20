# TwoPointFluxFVM



## [`FVMPhysics`](@ref)
This is an abstract type
which from which a user
data type can be derived which describes the physical
data of a given problem.

## [`FVMGraph`](@ref)

This is a weighted graph which represents edges and nodes
of a finite volume scheme. There are currently two
constructors:

[`FVMGraph(X::Array{Float64,1})`](@ref) for one-dimensional
domains and [`FVMGraph(X::Array{Float64,1},Y::Array{Float64,1})`](@ref)
for two-dimensional domains.

## [`TwoPointFluxFVMSystem`](@ref)

From instances of  [`FVMGraph`](@ref) and [`FVMPhysics`](@ref),
a [`TwoPointFluxFVMSystem`](@ref) which contains all the necessary
data for the solution of the nonlinear system described by them.




