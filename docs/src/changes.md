# Changes

## v0.8
- Replaced VoronoiFVM grid module by  [ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl)
- Moved grid generation, modification, plotting  over to ExtendableGrids
- Necessary changes in codes using VoronoiFVM:
  - Replace `grid.coord` by `coord` obtained via `coord=coordinates(grid)` or  `coord=grid[Coordinates]` after importing ExtendableGrids
  - Replace `VoronoiFVM.plot` by `ExtendableGrids.plot`.
  - In the plot method, Plotter is now a keyword argument
  - `VoronoiFVM.Grid()` now returns a ExtendableGrids.ExtendableGrid, 
    in fact it is just an alias to [ExtendableGrids.simplexgrid](https://j-fu.github.io/ExtendableGrids.jl/stable/search/?q=simplexgrid)
  - For using any methods on grids like cellmask! one nees to use `ExtendableGrids`
  - Subgrids now are of the same type `ExtendableGrids`,  views are currently defined for vectors only.
  
## v0.7 Feb 28 2019
- API modification:
  - __Breaking__:
    - `data` parameter passed to physics callbacks only if `Physics` object is created with `data` parameter.

      This makes the API more consistent in the case that parameters are just taken from the closure (the scope where the physics functions are defined)
      and no data `object` has been created.
    - Replace `node.coord[i]` by  `node[i]`.
    - Replace `edge.coordK[i]` by  `edge[i,1]`.
    - Replace `edge.coordL[i]` by  `edge[i,2]`.

      This now directly accesses the coordinate field of the grid and avoids copying of the coordinates
  - Backward compatible:  
    - No need for `viewK and viewL` in edge callbacks (they also make trouble with allocations...)
        - Replace `uk[i]` by `u[i,1]`
        - Replace `ul[i]` by `u[i,2]`
     - Replace `VoronoiFVM.DenseSystem(...)`  by `VoronoiFVM.System(..., unknown_storage=:dense)`
     - Replace `VoronoiFVM.SparseSystem(...)` by `VoronoiFVM.System(..., unknown_storage=:sparse)`

- No allocations anymore in assembly loop:
  - Replaced `ElasticArray` in `Grid` by normal one - this was the largest regression
  - Return `nothing` from mutating methods to avoid some allocations 
  - Indexing in `formfactors.jl` with `Int`

## v0.6.5 Jan 25 2019
- use updateindex! for matrix, depend on ExtendableSparse 0.2.6

## v0.6.4 2020-01-20
- Rearranged + commented boundary assembly loop
- Reworked + renamed some examples
- Document that unknowns doesn't initialize values.

## v0.6.3 2019-12-21 
remove xcolptrs call
Update dependency on ExtendableSparse

## v0.6.2 2019-12-20
Updated dependency list (Triangulate ^0.4.0)

## v0.6.1, 2019-12-17
- return "plotted" for being able to  place colormap
- require Triangulate >= 0.3.0

## v0.6.0, Dec 15 2019
- Removed Triangle submodule, depend on new Triangulate.jl Triangle wrapper
- link to source code in examples
- boundary_dirichlet! etc methods for setting boundary conditions
## v0.5.6 Dec 5 2019
- Bug fixes
- check triangle input for min 3 points
- check triangle edgelist for C_NULL
- voronoi plot
## v0.5.5 Dec 4 2019
- (Temporary) Copy of TriangleRaw as Triangle submodule. To be replaced by dependency on evisioned package
## v0.5.4 Dec 3 2019
- Re-enabled ElasticArrays in grid structure (for the time being)
- Added potkink example: this adds an inner boudary
## v0.5.3 Dec 1 2019
- triangle in optional submodule
- Modified API for plotting
   - Removed formal dependency on Plots and PyPlot
   - Use Plotter module as first parameter to plot methods  - replaces fvmplot
     and fvmpyplot functions. Use `VoronoiFVM.plot(PyPlot,...)` resp.  `VoronoiFVM.plot(Plots,...)`
   - No more complaints when package is used in environment with plots or pyplot installed
- Modified API for impedance

## v0.5.2 Nov 19, 2019
- Reorganized grid stuff
- Included triangle (after Ideas from TriangleMesh.jl)

## v0.5.1 Nov 13, 2019
- Fixed performance regression: AbstractArrays for Grid components were slow.
- Added handling of cylindrical coordinates

## V0.5, November 10, 2019
- Velocity projections
- Added edge handling to grid struct

## V0.4.2, November 6, 2019
- Replaced PyPlot by Plots
- Better and more examples

## V0.4, July 12, 2019
- Registered with Julia ecosystem
- Enhance Newton solver by embedding, exception handling
- Replace SparseMatrixCSC with ExtendableSparseMatrix
- fixed allocation issues in assembly
- assured that users get allocation stuff right via
  typed functions in physics structure
- more julianic API

## V0.3, April 9 2019
- Renamed from TwoPointFluxFVM to  VoronoiFVM
- Complete rewrite of assembly allowing sparse or dense matrix 
  to store degree of freedom information
    - Solution is a nnodes x nspecies sparse or dense matrix
    - The wonderful array interface of Julia still provides slicing
      etc in oder to access  species without need to write
      any bulk_solution stuff or whatever when using the sparse variant
- Re-export value() for debugging in physics functions
- Test function handling for flux calculation
- First working steps to impedance handling
- Abolished Graph in favor of  Grid, Graph was premature optimization...

## V0.2, Feb 20, 2019

- Changed signature of all callback functions:
  This also allows to pass user defined arrays etc. to the callback functions.
  In particular, velocity vectors can be passed this way.

  - Besides of `flux!()`, they now all have `node::VoronoiFVM.Node`
    as a second argument.

  - `flux!()` has `edge::VoronoiFVM.Edge` as a second argument

  - the `x` argument in `source!()` is omitted, the same data
     are now found in `node.coord`


- New method `edgelength(edge::VoronoiFVM.Edge)`
  
## V0.1, Dec. 2018

- Initial release
