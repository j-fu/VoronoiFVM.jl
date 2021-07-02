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


include("vfvm_xgrid.jl")
export cartesian!, circular_symmetric!, spherical_symmmetric!
export coordinates

include("vfvm_system.jl")
export unknowns
export enable_species!
export enable_boundary_species!
export update_grid!
export boundary_dirichlet!
export boundary_neumann!
export boundary_robin!
export value
export check_allocs!

include("vfvm_geometryitems.jl")
include("vfvm_formfactors.jl")
export meas,project
export unknown_indices
export edgevelocities,bfacevelocities


include("vfvm_newtoncontrol.jl")
export fixed_timesteps!,NewtonControl
export edgelength
export viewK,viewL,data


include("vfvm_solver.jl")
include("vfvm_diffeq_interface.jl")
export evolve!
export embed!
export solve!,solve



include("vfvm_integrals.jl")
include("vfvm_testfunctions.jl")
export integrate
export testfunction

include("vfvm_impedance.jl")
export impedance,freqdomain_impedance
export measurement_derivative


end

