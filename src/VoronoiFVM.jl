module VoronoiFVM

# Packages for meshing
using Triangulate

# Packages for Autodiff magic
using ForwardDiff, DiffResults
using IterativeSolvers
using DocStringExtensions
using VersionParsing


# These are in the standard distro
using SparseArrays
using ElasticArrays
using ExtendableSparse
using LinearAlgebra


using Printf


include("vfvm_physics.jl")
include("vfvm_grid.jl")
include("vfvm_functions.jl")
include("vfvm_newtoncontrol.jl")

include("vfvm_abstractsystem.jl")
include("vfvm_densesystem.jl")
include("vfvm_sparsesystem.jl")
include("vfvm_geometryitems.jl")
include("vfvm_subgridview.jl")

include("vfvm_solver.jl")
include("vfvm_oldapi.jl")
include("vfvm_testfunctions.jl")
include("vfvm_impedance.jl")
include("vfvm_plot.jl")



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
export num_nodes
export num_edges
export num_bfaces
export num_cells
export num_bfaceregions
export num_cellregions
export num_dof
export physics_data
export dim_space
export cellnode
export tridata
export prepare_edges!
export enable_species!
export enable_boundary_species!
export boundary_dirichlet!
export boundary_neumann!
export boundary_robin!
export cellmask!
export bfacemask!
export fvmplot
export fvmplot!
export bfacenode,nodecoord
export dof
export getdof
export setdof!
export value
export solve!
export embed!
export subgrid
export freqdomain_impedance
export measurement_derivative
export testfunction
export meas


export cartesian!, circular_symmetric!, spherical_symmmetric!
export edgevelocities
export geomspace
export TokenStream,gettoken, expecttoken,trytoken


# deprecated
export edgelength
export viewK,viewL


end

