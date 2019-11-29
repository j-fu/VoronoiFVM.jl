module VoronoiFVM

import  Pkg

# Packages for Autodiff magic
using ForwardDiff, DiffResults
using IterativeSolvers
using DocStringExtensions


# These are in the standard distro
using SparseArrays
using ElasticArrays
using ExtendableSparse
using LinearAlgebra


using Printf

installed(pkg)=haskey(Pkg.installed(),pkg)

if !installed("TriangleRaw")
    include("triangle/triangle.jl")
end


include("vfvm_physics.jl")
include("vfvm_grid.jl")
include("vfvm_functions.jl")
include("vfvm_newtoncontrol.jl")
include("vfvm_system.jl")
include("vfvm_solver.jl")
include("vfvm_testfunctions.jl")
include("vfvm_impedance.jl")



if installed("Plots")
    include("vfvm_plots.jl")
end



if installed("PyPlot")
    include("vfvm_pyplot.jl")
end

export installed


export unknowns
export fbernoulli
export fbernoulli_pm
export integrate



export glue
export num_nodes
export num_edges
export num_bfaces
export num_cells
export num_bfaceregions
export num_cellregions
export num_dof
export tridata
export prepare_edges!
export enable_species!
export enable_boundary_species!
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
export data
export edgelength
export viewK,viewL
export cartesian!, circular_symmetric!, spherical_symmmetric!
export edgevelocities
export geomspace
export TokenStream,gettoken, expecttoken,trytoken
export TriangulateIO, triangulate
end

