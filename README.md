VoronoiFVM.jl
===============

[![Build Status](https://img.shields.io/travis/j-fu/VoronoiFVM.jl/master.svg?label=Linux+MacOSX)](https://travis-ci.org/j-fu/VoronoiFVM.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://j-fu.github.io/VoronoiFVM.jl/dev)


Author: J. Fuhrmann


Solver for coupled nonlinear partial differential equations based on the Voronoi finite volume method.


This Julia package merges the ideas behind [pdelib](http://www.wias-berlin.de/software/pdelib/?lang=0)/fvsys with the multiple dispatch paradigm of Julia. It allows to avoid the implementation of function derivatives.  It instead uses [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and [DiffResults](https://github.com/JuliaDiff/DiffResults.jl) to evaluate user functions along with their jacobians.


It requires Julia 1.x

The package contains a submodule Triangle which is a current copy of the provisional package [TriangleRaw.jl](https://github.com/j-fu/TriangleRaw)
which attempts to consolidate the use of triangle based on ideas from [TriangleMesh.jl](https://github.com/konsim83/TriangleMesh.jl)
and [Triangle.jl](https://github.com/cvdlab/Triangle.jl).

