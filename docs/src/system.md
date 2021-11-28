# System

The computational grid is required is assumed to correspond to a domain
``\Omega=\cup_{r=1}^{n_\Omega} \Omega_r`` 

with boundary  ``\partial\Omega=\Gamma=\cup_{b=1}^{n_\Gamma} \Gamma_b``.

The subdomains ``\Omega_r`` are called "regions" and the boundary
subdomains ``\Gamma_b`` are called "boundary regions".

On this complex of domains "lives"  a number of species which are either
attached to a number of regions or to a number of boundary regions.

Grids for VoronoiFVM are managed by the package 
[ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl).


```@autodocs
Modules = [VoronoiFVM]
Pages = ["vfvm_system.jl"]
```

