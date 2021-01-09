VoronoiFVM.jl
===============

[![Build status](https://github.com/j-fu/VoronoiFVM.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/j-fu/VoronoiFVM.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://j-fu.github.io/VoronoiFVM.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://j-fu.github.io/VoronoiFVM.jl/dev)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3529808.svg)](https://doi.org/10.5281/zenodo.3529808)


Solver for coupled nonlinear partial differential equations based on the Voronoi finite volume method.


This Julia package merges the ideas behind [pdelib](http://www.wias-berlin.de/software/pdelib/?lang=0)/fvsys with the multiple dispatch paradigm of Julia. It allows to avoid the implementation of function derivatives.  It instead uses [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and [DiffResults](https://github.com/JuliaDiff/DiffResults.jl) to evaluate user functions along with their jacobians.


It requires Julia 1.x

API changes:
- [v0.10](https://j-fu.github.io/VoronoiFVM.jl/v0.10/changes)
- [v0.9](https://j-fu.github.io/VoronoiFVM.jl/v0.9/changes/) (some are breaking)
- [v0.8](https://j-fu.github.io/VoronoiFVM.jl/v0.8/changes/) (some are breaking)
- [v0.7](https://j-fu.github.io/VoronoiFVM.jl/v0.7/changes/) (some are breaking)


