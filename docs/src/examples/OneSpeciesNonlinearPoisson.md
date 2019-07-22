# 1D Nonlinear Poisson equation with one species

Solve the nonlinear Poisson equation

```math
-\nabla 0.01 \nabla u^2 + u^2 = 0.0001x
```
in $\Omega=(0,1)$ with boundary condition $u(0)=1$ and $u(1)=0.5$.

````@eval
using Markdown
Markdown.parse("""
```julia
$(read("../../../examples/OneSpeciesNonlinearPoisson.jl",String))
```
""")
````
