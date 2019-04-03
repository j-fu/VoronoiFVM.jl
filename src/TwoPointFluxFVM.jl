module TwoPointFluxFVM

include("fvmgrid.jl")
include("functions.jl")
include("fvmnewtoncontrol.jl")
include("tools.jl")
include("fvmpyplot.jl")
include("twopointfluxfvmsystem.jl")


export unknowns
export bulk_unknowns
export boundary_unknowns
export fbernoulli
export fbernoulli_pm
export solve
export integrate




export nnodes
export nbfaces
export ncells
export add_species
export add_boundary_species
export copy
export cellmask!
export fvmplot
export dof
export getdof
export setdof!
export value

end

