VoronoiFVM.jl
===============

[![Build status](https://github.com/j-fu/VoronoiFVM.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/j-fu/VoronoiFVM.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://j-fu.github.io/VoronoiFVM.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://j-fu.github.io/VoronoiFVM.jl/dev)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3529808.svg)](https://doi.org/10.5281/zenodo.3529808)
[![Zulip Chat](https://img.shields.io/badge/zulip-join_chat-brightgreen.svg)](https://julialang.zulipchat.com/#narrow/stream/379007-voronoifvm.2Ejl)


Solver for coupled nonlinear partial differential equations (elliptic-parabolic conservation laws) based on the Voronoi finite volume method.
It uses automatic differentiation via [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) and [DiffResults.jl](https://github.com/JuliaDiff/DiffResults.jl) to evaluate user functions along with their jacobians and calculate derivatives of solutions with respect to their parameters.

## Recent changes
Please look up the list of recent [changes](https://j-fu.github.io/VoronoiFVM.jl/stable/changes)

## Accompanying packages
- [ExtendableSparse.jl](https://github.com/j-fu/ExtendableSparse.jl): convenient and efficient sparse matrix assembly
- [ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl): unstructured grid management library
- [SimplexGridFactory.jl](https://github.com/j-fu/SimplexGridFactory.jl): unified high level  mesh generator interface
- [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl):  Julia wrapper for the [Triangle](https://www.cs.cmu.edu/~quake/triangle.html) triangle mesh generator by J. Shewchuk
- [TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl):  Julia wrapper for the [TetGen](http://www.tetgen.org) tetrahedral mesh generator by H. Si.
- [GridVisualize.jl](https://github.com/j-fu/GridVisualize.jl): grid and function visualization related to ExtendableGrids.jl
- [PlutoVista.jl](https://github.com/j-fu/PlutoVista.jl): backend for [GridVisualize.jl](https://github.com/j-fu/GridVisualize.jl) for use in Pluto notebooks.

VoronoiFVM.jl and most of these packages are  part of the meta package [PDELib.jl](https://github.com/WIAS-BERLIN/PDELib.jl).



## Some alternatives
- [ExtendableFEM.jl](https://github.com/chmerdon/ExtendableFEM.jl): finite element library implementing gradient robust FEM
  from the same package base by Ch. Merdon
- [SkeelBerzins.jl](https://github.com/gregoirepourtier/SkeelBerzins.jl): a Julian variation on Matlab's `pdepe` API
- [Trixi.jl](https://github.com/trixi-framework/Trixi.jl):  numerical simulation framework for hyperbolic conservation laws 
- [GridAP.jl](https://github.com/gridap/Gridap.jl) Grid-based approximation of partial differential equations in Julia
- [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) Finite element toolbox for Julia
- [FinEtools.jl](https://github.com/PetrKryslUCSD/FinEtools.jl)  Finite element tools for Julia
- [FiniteVolumeMethod.jl](https://github.com/DanielVandH/FiniteVolumeMethod.jl/) Finite volumes with [Donald boxes](https://sciml.github.io/FiniteVolumeMethod.jl/dev/math/#Control-volumes)

## Some projects and packages using VoronoiFVM.jl

- [RfbScFVM: Performance prediction of flow battery vells](https://github.com/Isomorph-Electrochemical-Cells/RfbScFVM)
- [ChargeTransport.jl: Drift diffusion simulator for semiconductor devices](https://github.com/PatricioFarrell/ChargeTransport.jl)
- [MosLab.jl: From semiconductor to transistor level modeling in Julia](https://github.com/Rapos0/MOSLab.jl)
- [LiquidElectrolytes.jl: Generalized Nernst-Planck-Poisson model for liquid electrolytes](https://github.com/j-fu/LiquidElectrolytes.jl)


## Citation

If you use this package in your work, please cite it according to [CITATION.cff](https://raw.githubusercontent.com/j-fu/VoronoiFVM.jl/master/CITATION.cff)
