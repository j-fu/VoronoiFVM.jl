# Postprocessing

## Plotting

Plotting can be performed using the package [GridVisualize.jl](https://github.com/j-fu/GridVisualize.jl).
This package extends the API with a couple of methods:
    
```@docs
GridVisualize.gridplot
GridVisualize.gridplot!
GridVisualize.scalarplot
GridVisualize.scalarplot!
```

## Norms & volumes
```@docs
LinearAlgebra.norm
lpnorm
l2norm
w1pseminorm
h1seminorm
w1pnorm
h1norm
lpw1pseminorm
l2h1seminorm
lpw1pnorm
l2h1norm
nodevolumes
```

## Solution integrals
```@docs
integrate(system::VoronoiFVM.AbstractSystem, F::Function, U::AbstractMatrix; boundary = false)
VoronoiFVM.edgeintegrate
```

## Nodal flux reconstruction
```@docs
nodeflux
```

## Boundary flux calculation
```@autodocs
Modules = [VoronoiFVM]
Pages = ["vfvm_testfunctions.jl"]
```

## Impedance calculatiom
Impedance calculation can be seen as a postprocessing step
after the solution of the unexcited stationary system.


```@autodocs
Modules = [VoronoiFVM]
Pages = ["vfvm_impedance.jl"]
```

