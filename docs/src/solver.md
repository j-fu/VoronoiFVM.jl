# Solution methods

## Solve method

```@docs
solve(system::VoronoiFVM.AbstractSystem; kwargs...)
solve(inival, system::VoronoiFVM.AbstractSystem, times; kwargs...)
solve(inival, system::VoronoiFVM.AbstractSystem; kwargs...)
solve!(solution,inival, system::VoronoiFVM.AbstractSystem; kwargs...)
evolve!(solution, inival, system::VoronoiFVM.AbstractSystem, times; kwargs...)
embed!(solution, inival,system::VoronoiFVM.AbstractSystem; kwargs...)
``` 


## Solver control
```@docs 
SolverControl
NewtonControl
```
    

## [Interface to DifferentialEquations.jl](@id diffeq)

```@autodocs
Modules = [VoronoiFVM]
Pages = ["vfvm_diffeq_interface.jl"]
```

