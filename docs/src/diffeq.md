# Interface to DifferentialEquations.jl

Since version 0.9 VoronoiFVM provides an experimental interface to the ODE solvers
in DifferentialEquations.jl. The general pattern is as follows:

```julia
using VoronoiFVM
using DifferentialEquations

sys=VoronoiFVM.System(...)

tspan = (t0,tend)
    
ode_func = ODEFunction(eval_rhs!,
                       jac=eval_jacobian!,
                       jac_prototype=jac_prototype(sys),
                       mass_matrix=mass_matrix(sys))
    
problem = ODEProblem(ode_func,vec(inival),tspan,sys)
    sol = DifferentialEquations.solve(problem,Rodas5())
```


Some caveats with this evolving implementation:

- Ensure that you handle problems wit  "trivial" storage function:
```julia
    storage!(f,u,node)= f.=u
```

- Eventually it is possible to generalize this to multiplication of `u` by species and/or
  region dependent constants, but this has to be implemented  with the mass matrix


- The choice of ODE methods is limited to those [able to handle mass matrices](https://diffeq.sciml.ai/stable/solvers/dae_solve/#OrdinaryDiffEq.jl-(Mass-Matrix))



```@autodocs
Modules = [VoronoiFVM]
Pages = ["vfvm_diffeq_interface.jl"]
```

