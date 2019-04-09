module TwoPointFluxFVM

include("fvmgrid.jl")
include("functions.jl")
include("fvmnewtoncontrol.jl")
include("tools.jl")
include("fvmpyplot.jl")
include("twopointfluxfvmsystem.jl")
include("testfunctions.jl")
include("fvmimpedance.jl")


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

