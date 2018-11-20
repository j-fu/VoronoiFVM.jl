TwoPointFluxFVM
===============

Author: J. Fuhrmann


Solver for coupled nonlinear PDEs based on two point flux finite volume method on admissible grids.


This Julia package merges the ideas behind pdelib/fvsys with the multiple dispatch paradigm of Julia. It allows to avoid the implementation of function derivatives. It instead uses the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and [DiffResult](https://github.com/JuliaDiff/DiffResult.jl) and to evaluate user fucntions along with their jacobians.


So far this is merely a design study to learn and evaluate Julia.  It is however aimed to be feasible at least for small projects.

It requires Julia 1.0.

## Documentation

It currently resides at [WIAS Berlin](https://www.wias-berlin.de/people/fuhrmann/TwoPointFluxFVM/)

## Installation 

So far, the package has not been registered with Julia.

Steps:

   1. Create a directory for Julia packages which have not been registered, let us denote this
      as PKG_DIR
    
   2. cd to PKG_DIR and clone this repository:

      Either from WIAS rhodecode server
````
     git clone  https://repos.wias-berlin.de/users/fuhrmann/projects/julia-packages/TwoPointFluxFVM
````
     Or from github:
````
     git clone  https://github.com/j-fu/TwoPointFluxFVM.jl TwoPointFluxFVM
````
   
   3. Add the following line to  the file `.julia/config/startup.jl` in your  home directory (create this file it does not exist). PKG_DIR must be the full path name of the directory.
````
     push!(LOAD_PATH, "PKG_DIR")
````

Now, `import TwoPointFluxFVM` should work in Julia scripts.
