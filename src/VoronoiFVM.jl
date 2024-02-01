"""
$(EXPORTS)
"""
module VoronoiFVM
using Printf
using DocStringExtensions
using LinearAlgebra
using SparseArrays
import CommonSolve, SciMLBase


using InteractiveUtils
using BandedMatrices
# using MultidiagonalMatrices

using LinearSolve
using Statistics

using ForwardDiff
using DiffResults
using JLD2
using RecursiveArrayTools

using ExtendableSparse
using ExtendableGrids

using StaticArrays
using SparseDiffTools
using Symbolics
using Random
using RecursiveFactorization

include("vfvm_physics.jl")
include("vfvm_functions.jl")
export fbernoulli
export fbernoulli_pm
export inplace_linsolve!

include("vfvm_densesolution.jl")
include("vfvm_sparsesolution.jl")
export num_dof
export dof
export getdof
export setdof!

include("vfvm_history.jl")
export NewtonSolverHistory, TransientSolverHistory, details

include("vfvm_transientsolution.jl")
export transient_solution, TransientSolution


include("vfvm_xgrid.jl")
export cartesian!, circular_symmetric!, spherical_symmetric!
export coordinates

"""
$(TYPEDEF)
    
Abstract type for finite volume system structure.
"""
abstract type AbstractSystem{Tv <: Number, Tc <: Number, Ti <: Integer, Tm <: Integer} end
include("vfvm_geometryitems.jl")
include("vfvm_assemblydata.jl")
include("vfvm_system.jl")
export unknowns
export num_species
export partitioning, Equationwise
export enable_species!
export enable_boundary_species!
export enable_discontinuous_species!
export update_grid!
export boundary_dirichlet!
export boundary_neumann!
export boundary_robin!
export ramp
export value
export physics!
export history, history_summary, history_details
export evaluate_residual_and_jacobian
export edgelength
export viewK, viewL, data
export hasoutflownode,isoutflownode, outflownode

include("vfvm_formfactors.jl")
export meas, project
export unknown_indices
export edgevelocities, bfacevelocities, bfacenodefactors
export time, region, embedparam, parameters
export calc_divergences

include("vfvm_solvercontrol.jl")
export fixed_timesteps!, NewtonControl, SolverControl
include("vfvm_linsolve.jl")
export DirectSolver, GMRESIteration, CGIteration, BICGstabIteration, NoBlock, EquationBlock, PointBlock


include("vfvm_exceptions.jl")
include("vfvm_assembly.jl")
include("vfvm_solver.jl")
export evolve!
export embed!
export solve!, solve

include("vfvm_postprocess.jl")
export nodeflux
export integrate
export l2norm,lpnorm
export w1pseminorm,h1seminorm
export w1pnorm,h1norm
export lpw1pnorm,l2h1norm
export lpw1pseminorm,l2h1seminorm
export nodevolumes

include("vfvm_testfunctions.jl")
export testfunction
export TestFunctionFactory

include("vfvm_quantities.jl")
export ContinuousQuantity
export DiscontinuousQuantity
export InterfaceQuantity
export subgrids, views

include("vfvm_impedance.jl")
export impedance, freqdomain_impedance
export measurement_derivative

include("vfvm_diffeq_interface.jl")
export eval_rhs!, eval_jacobian!, mass_matrix, prepare_diffeq!

include("gridvisualize.jl")
export plothistory
#include("precompile.jl")

end
