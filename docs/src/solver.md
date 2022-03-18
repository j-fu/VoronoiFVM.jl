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
For transient problems, as an alternative to the use of the built-in implicit Euler method, (stiff) ODE solvers from 
[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)  can be used.

The dependency of the code on  [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) is optional.
It is handeled via [Requires.jl](https://github.com/JuliaPackaging/Requires.jl). This means that it becomes available as
as soon as [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) is made available via `Pkg.add()`
and used in the code.

The interface just provides two methods: creation of an [ODEProblem](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems) from a [`VoronoiFVM.System`](@ref) and a reshape method
which turns the output of the ode solver into a [`TransientSolution`](@ref).

The basic usage pattern is as follows: replace the call to the built-in solver
```julia
sol=VoronoiFVM.solve(sys; times=(t0,t1), inival=inival)
```
by
```julia
problem = ODEProblem(sys,inival,(t0,tend))
odesol = DifferentialEquations.solve(problem)
sol=reshape(odesol,sys)
```


```@docs
DifferentialEquations.ODEProblem(::VoronoiFVM.AbstractSystem, inival, tspan, callback)
Base.reshape(::DifferentialEquations.AbstractDiffEqArray,::VoronoiFVM.AbstractSystem)
```


## Built-in solver control
```@docs 
SolverControl
NewtonControl
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
