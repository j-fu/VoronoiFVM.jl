"""
$(README)

$(EXPORTS)
"""
module VoronoiFVM



# Packages for Autodiff magic
using ForwardDiff, DiffResults
using IterativeSolvers
using DocStringExtensions


# These are in the standard distro
using SparseArrays
using ExtendableSparse
using LinearAlgebra

using SparseDiffTools
using SparseArrays
using SparsityDetection

using Printf

# Packages for meshing -> remove from here!
# using Triangulate

using ExtendableGrids
include("vfvm_xgrid.jl")
#include("vfvm_grid.jl")



include("vfvm_physics.jl")
include("vfvm_functions.jl")
include("vfvm_newtoncontrol.jl")

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


export isplots
export ispyplot
export ispyplotter

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
export value
export solve!
export evolve!
export embed!
export freqdomain_impedance
export measurement_derivative
export testfunction
export meas
export unknown_indices

export cartesian!, circular_symmetric!, spherical_symmmetric!
export edgevelocities

# deprecated
export edgelength
export viewK,viewL,data

end

