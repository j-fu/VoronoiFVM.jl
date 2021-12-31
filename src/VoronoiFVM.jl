"""
$(README)

$(EXPORTS)
"""
module VoronoiFVM

using Printf
using DocStringExtensions
using LinearAlgebra
using SparseArrays
using SuiteSparse

using Parameters
using Statistics

using ForwardDiff
using DiffResults
using IterativeSolvers
using JLD2
using StaticArrays
using SparseDiffTools
using SparsityDetection
using RecursiveArrayTools

using ExtendableSparse
using ExtendableGrids


include("vfvm_physics.jl")
include("vfvm_functions.jl")
export fbernoulli
export fbernoulli_pm


include("vfvm_densesolution.jl")
include("vfvm_sparsesolution.jl")
export num_dof
export dof
export getdof
export setdof!

include("vfvm_transientsolution.jl")
export transient_solution,TransientSolution

include("vfvm_history.jl")
export NewtonSolverHistory, TransientSolverHistory, details

include("vfvm_xgrid.jl")
export cartesian!, circular_symmetric!, spherical_symmmetric!
export coordinates

include("vfvm_system.jl")
export unknowns
export enable_species!
export enable_boundary_species!
export enable_discontinuous_species!
export update_grid!
export boundary_dirichlet!
export boundary_neumann!
export boundary_robin!
export ramp
export value
export check_allocs!
export physics!
export history,history_summary,history_details

include("vfvm_geometryitems.jl")
include("vfvm_formfactors.jl")
export meas,project
export unknown_indices
export edgevelocities,bfacevelocities
export time,region,embedparam

include("vfvm_solvercontrol.jl")
export fixed_timesteps!,NewtonControl,SolverControl
export edgelength
export viewK,viewL,data


include("vfvm_solver.jl")
include("vfvm_diffeq_interface.jl")
export evolve!
export embed!
export solve!,solve


include("vfvm_postprocess.jl")
export nodeflux
include("vfvm_testfunctions.jl")
export integrate
export testfunction

include("vfvm_quantities.jl")
export ContinuousQuantity
export DiscontinuousQuantity
export InterfaceQuantity


include("vfvm_impedance.jl")
export impedance,freqdomain_impedance
export measurement_derivative

include("gridvisualize.jl")

include("precompile.jl")


end

