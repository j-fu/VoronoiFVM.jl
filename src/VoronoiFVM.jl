module VoronoiFVM

include("vfvm_grid.jl")
include("vfvm_functions.jl")
include("vfvm_newtoncontrol.jl")
include("vfvm_tools.jl")
include("vfvm_pyplot.jl")
include("vfvm_system.jl")
include("vfvm_testfunctions.jl")
include("vfvm_impedance.jl")


export unknowns
export fbernoulli
export fbernoulli_pm
export integrate




export num_nodes
export num_bfaces
export num_cells
export add_species
export add_boundary_species
export cellmask!
export fvmplot
export dof
export getdof
export setdof!
export value
export solve
export subgrid
export testfunction


end

