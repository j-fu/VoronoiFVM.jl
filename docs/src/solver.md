# Solution methods

## Solve method

```@docs
solve(system::VoronoiFVM.AbstractSystem; kwargs...)
solve(inival, system::VoronoiFVM.AbstractSystem, times; kwargs...)
solve(inival, system::VoronoiFVM.AbstractSystem; kwargs...)
solve!(solution,inival, system::VoronoiFVM.AbstractSystem; kwargs...)
``` 
## Solver history
If `log` is set to true in `solve`, the history of newton iterations and  time/embedding
steps is recorded and. For the respective previous solution step it can be obtained via
`history(system)`.

```@docs
NewtonSolverHistory
TransientSolverHistory
Base.summary
details
history
history_details
history_summary
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

