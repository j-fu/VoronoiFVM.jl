# Plotting
In order to avoid drawing in of heavy dependencies for VoronoiFVM,
the plot methods for grids and grid functions defined in this package
have as their first argument the module of the plotting package used.

Currently, PyPlot and Plots are supported. Similar schemes in the future
might work for Makie and other packages.


```@autodocs
Modules = [VoronoiFVM]
Pages = ["vfvm_plot.jl"]
Order = [:function]
```
