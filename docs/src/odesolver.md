# [OrdinaryDiffEq.jl solver](@id diffeq)


For transient problems, as an alternative to the use of the built-in implicit Euler method, (stiff) ODE solvers from 
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl)  can be used.

The interface just provides two methods: creation of an [ODEProblem](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems) from a [`VoronoiFVM.System`](@ref) and a reshape method
which turns the output of the ode solver into a [`TransientSolution`](@ref).

The basic usage pattern is as follows: use [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) and replace the call to the built-in solver
```julia
sol=solve(sys; times=(t0,t1), inival=inival)
```
by
```julia
problem = ODEProblem(sys,inival,(t0,tend))
odesol = solve(problem, solver)
sol=reshape(odesol,sys)
```
Here, `solver` is one of the aforementioned ODE solvers.  If the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
package is used, the `solver` parameter is omitted, and some default is chosen.

```@docs
SciMLBase.ODEFunction(sys::VoronoiFVM.AbstractSystem; jacval=unknowns(sys,0), tjac=0)
SciMLBase.ODEProblem(::VoronoiFVM.AbstractSystem, inival, tspan, callback)
Base.reshape(::RecursiveArrayTools.AbstractDiffEqArray,::VoronoiFVM.AbstractSystem)
```


## Choice of ODE/DAE solver

As this package interfaces to  the PDE solver package [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl),
a general advise is to choose methods able to handle [stiff problems](https://diffeq.sciml.ai/stable/solvers/ode_solve/#Stiff-Problems).
Moreover, often, discretized PDE systems (e.g. containing elliptic equations) are differential agebraic equation (DAE) systems 
which should be solved by [DAE solvers](https://diffeq.sciml.ai/stable/solvers/dae_solve/).


