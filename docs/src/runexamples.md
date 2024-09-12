About the examples
==================
The examples have been designed with the following issues in mind:
- they run from the Julia REPL
- each example is a Julia module named similar to the basename of the example file.
- an example can be used as the starting point for a project 
- the examples at the same time comprise the test suite for VoronoiFVM.

Since the creation of these examples, the API has been updated and simplified.


## Running the examples
Plotting is performed using the GridVisualize.jl package which interfaces PyPlot.jl,
Plots.jl, Makie.jl.

In order to run `ExampleXXX`, perform the following steps:

- Download the example file (e.g. via the source code link at the top)
- Call Julia with  an Julia environment which contains `VoronoiFVM.jl`, `ExtendableGrids.jl`, `GridVisualize.jl` 
  and e.g. `PyPlot.jl`
- `include("ExampleXXX.jl")`
- Run the example via `ExampleXXX.main(Plotter=PyPlot)`

Due to the encapsulation into modules, you can load as many examples as you like.

If you want to modify the example, consider using `Revise.jl` and `includet`. 


## Performance with closures

VoronoiFVM provides two flavors of calbacks for constitutive
functions: 

- Callbacks with `data` parameoter. `data` is declared as part of `Physics` and
  passed down to the callbacks
- Callbacks without `data` parameter. Here, the parameters of the 
  physics callbacks are accessed via closures, i.e. from within the 
  scope of the definition of the particular function.

While the  second method is  very convenient to  use, it comes  with a
serious performance pitfall: if a  variable in the closure is assigned
twice, Julia becomes unsure about  it's type and therefore "boxes" it,
i.e. it  creates a wrapper struct  around the variable value  which is
able to track  its potentially changing type.  The serious consequence
of this is that assignments to a boxed variable lead to allocations, which
are a serious performance hit if they occur in loops over grid nodes
or edges.

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
