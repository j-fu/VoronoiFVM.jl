TwoPointFluxFVM
===============

Author: J. Fuhrmann


Solver for coupled nonlinear PDEs based on two point flux finite volume
method on admissible grids.


This Julia package merges the ideas behind pdelib/fvsys with
the multiple dispatch paradigm of Julia. It allows to avoid
the implementation of function derivatives. It instead uses the
[ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and 
[DiffResult](https://github.com/JuliaDiff/DiffResult.jl) and 
to evaluate user fucntions along with their jacobians.


So far this is merely a design study to learn and evaluate Julia.
It is however aimed to be feasible at least for small projects.



