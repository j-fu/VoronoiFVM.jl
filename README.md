TwoPointFluxFVMSystem
=====================

Solver for coupled nonlinear PDEs based on two point flux finite volume
method on admissible grids.


This Julia package merges the ideas behind pdelib/fvsys with
the multiple dispatch paradigm of Julia. It allows to avoid
the implementation of function derivatives. It instead uses the
ForwardDiff and DiffResult packages in order to calculate
jacobians of functions.


