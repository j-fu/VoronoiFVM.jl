module TwoPointFluxFVM

include("fvmgrid.jl")
include("functions.jl")
include("fvmnewtoncontrol.jl")
include("tools.jl")
include("fvmpyplot.jl")
include("twopointfluxfvmsystem.jl")


export unknowns
export fbernoulli
export fbernoulli_pm
export solve
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
export FVMSubGrid
export FVMGrid
export FVMNewtonControl
export SparseFVMSystem
export solve

end

