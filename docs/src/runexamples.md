About the examples
==================
The examples have been designed with the following issues in mind:
- they run from the Julia REPL
- each example is a Julia module named similar to the basename of the example file.
- an example can be used as the starting point for a project 
- the examples at the same time comprise the test suite for VoronoiFVM.


## Running the examples
Plotting is performed using the GridVisualize.jl package which interfaces PyPlot.jl,
Plots.jl, Makie.jl.

In order to run `ExampleXXX`, peform the following steps:

- Download the example file (e.g. via the source code link at the top)
- Call Julia with  an Julia environment which contains `VoronoiFVM.jl`, `ExtendableGrids.jl`, `GridVisualize.jl` 
  and e.g. `PyPlot.jl`
- `include("ExampleXXX.jl")`
- Run the example via `ExampleXXX.main(Plotter=PyPlot)`

Due to the encapsulation into modules, you can load as many examples as you like.

If you want to modifiy the example, consider using `Revise.jl` and `includet`. 


## Performance with  closures

VoronoiFVM provides two flavors of calbacks for constitutive
functions: one, where `data` is declared as part of `Physics` and
passed down to the callbacks, and one, where the parameters of the
physics callbacks are passed to them from the scope of the definition
of the functions, via so-called closures.

While the second method is very convenient to use, it comes with a
serious performance pitfall: if a variable is assigned twice, Julia
becomes unsure about it's type and therefore "boxes" it, i.e. it creates
a wrapper struct around the variable value which is able to track its
potentially changing type. The serious consequence of this is that 
assignments to a boxed variable lead to allocations.

This behaviour is explained in the [Julia documentation](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured).

Here is an example which comes close to the situation in VoronoiFVM:

```@example
function ttype_boxed(n)
    u=rand(n)
    v=similar(u)
    a=2.0
    a=3.0
    dostuff(u)=a*u
    @allocated map!(dostuff,v,u)
end
ttype_boxed(10) # hide
ttype_boxed(10)
```

The remedy is to type-annotate variables from closures:
```@example
function ttype_annotated(n)
    u=rand(n)
    v=similar(u)
    a::Float64=2.0
    a=3.0
    dostuff(u)=a*u
    @allocated map!(dostuff,v,u)
end
ttype_annotated(10) # hide
ttype_annotated(10)
```

VoronoiFVM provides a mechanism to check for this situation with the
[`check_allocs!`](@ref) method. After creating a system, set
```julia
check_allocs!(system,true)
```
and allocations inner assembly loops will throw an error.

To set this behavior as default for all systems you create, invoke 
```julia
ENV["VORONOIFVM_CHECK_ALLOCS"]="true"
```
before using `VoronoiFVM`, e.g. put it  into your `startup.jl`
