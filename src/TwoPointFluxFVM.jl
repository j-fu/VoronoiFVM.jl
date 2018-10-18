module TwoPointFluxFVM

include("fvmgraph.jl")
include("functions.jl")
include("fvmparameters.jl")
include("fvmnewtoncontrol.jl")
include("twopointfluxfvmsystem.jl")



export FVMGraph
export TwoPointFluxFVMSystem
export FVMParameters
export DefaultFVMParameters
export FVMNewtonControl
export unknowns
export bulk_unknowns
export boundary_unknowns
export fbernoulli
export fbernoulli_pm
export solve
export integrate
export Dirichlet
export @AddDefaultFVMParameters
end

