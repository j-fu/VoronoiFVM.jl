# 1D Nonlinear Poisson equation with two species
Solve the nonlinear Poisson equation

```math
-\nabla (0.01+2u_2)\nabla u_1 + u_1u_2= 0.0001(0.01+x)
```

```math
-\nabla (0.01+2u_1)\nabla u_2 -+ u_1u_2 = 0.0001(1.01-x)
```


in $\Omega=(0,1)$ with boundary condition $u_1(0)=1$, $u_2(0)=0$ and $u_1(1)=1$, $u_2(1)=1$.

````@eval
using Markdown

Markdown.parse("""
```julia
$(read("../../../examples/TwoSpeciesNonlinearPoisson.jl",String))
```
""")
````
