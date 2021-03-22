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


include("vfvm_xgrid.jl")
export cartesian!, circular_symmetric!, spherical_symmmetric!



include("vfvm_physics.jl")
include("vfvm_functions.jl")
include("vfvm_newtoncontrol.jl")
export fixed_timesteps!,NewtonControl

include("vfvm_densesolution.jl")
include("vfvm_sparsesolution.jl")
include("vfvm_transientsolution.jl")

include("vfvm_abstractsystem.jl")
include("vfvm_densesystem.jl")
include("vfvm_sparsesystem.jl")

include("vfvm_geometryitems.jl")
include("vfvm_formfactors.jl")



include("vfvm_solver.jl")
include("vfvm_testfunctions.jl")
include("vfvm_impedance.jl")



export unknowns
export fbernoulli
export fbernoulli_pm
export integrate


export FVMSystem
export FVMPhysics

export glue
export coordinates
export num_dof
export physics_data
export dim_space
export cellnode
export enable_species!
export enable_boundary_species!
export update_grid!
export boundary_dirichlet!
export boundary_neumann!
export boundary_robin!
export dof
export getdof
export setdof!
export value,partials,npartials
export evolve!
export embed!
export freqdomain_impedance
export measurement_derivative
export testfunction
export meas,project
export unknown_indices

export edgevelocities

export edgelength
export viewK,viewL,data

include("vfvm_diffeq_interface.jl")
export eval_rhs!,eval_jacobian!,mass_matrix,jac_prototype
export transient_solution,TransientSolution
export solve!,solve


end

