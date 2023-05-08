# Solution methods

## Built-in solver
The package comes with a built-in solve method which solves 
stationary problems, simple homotopy embedding problems and transient problems 
via the implicit Euler method.  This solver and its default parameters are tuned for robustness,
possibly at the expense of solution speed. Careful tuning of the parameters, or -- in the case of transient problems --
the choice of the [DifferentialEquations.jl solver](@ref diffeq) can significantly improve the performance.

```@docs
VoronoiFVM.solve(system::VoronoiFVM.AbstractSystem; kwargs...)
``` 

## [DifferentialEquations.jl solver](@id diffeq)
This solver has been being outsourced into the glue package VoronoiFVMDiffEq.jl.


## Built-in solver control
```@docs 
SolverControl
NewtonControl
```

## Linear Solver strategies
```@docs
VoronoiFVM.SolverStrategies
VoronoiFVM.SolverStrategies.direct_umfpack
VoronoiFVM.SolverStrategies.gmres_umfpack
VoronoiFVM.SolverStrategies.gmres_eqnblock_umfpack
VoronoiFVM.SolverStrategies.gmres_iluzero
VoronoiFVM.SolverStrategies.gmres_eqnblock_iluzero
VoronoiFVM.SolverStrategies.gmres_pointblock_iluzero
```


## Built-in solver history handling
If `log` is set to true in `solve`, the history of newton iterations and  time/embedding
steps is recorded and. For the respective previous solution step it can be obtained via
`history(system)`.

```@docs
NewtonSolverHistory
TransientSolverHistory
Base.summary(::NewtonSolverHistory)
Base.summary(::TransientSolverHistory)
details
history
history_details
history_summary
```

## Matrix extraction
For testing and teaching purposes, one can obtain residual and linearization at a given vector of unknowns

```@docs
evaluate_residual_and_jacobian
```

## Legacy API
During the development of the code, a number of API variants have been developed which are supported for backward compatibility.

```@docs
VoronoiFVM.solve(inival, system::VoronoiFVM.AbstractSystem, times; kwargs...)
VoronoiFVM.solve(inival, system::VoronoiFVM.AbstractSystem; kwargs...)
VoronoiFVM.solve!(solution,inival, system::VoronoiFVM.AbstractSystem; kwargs...)
``` 
