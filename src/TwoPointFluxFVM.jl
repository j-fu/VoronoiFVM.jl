module TwoPointFluxFVM

include("fvmgraph.jl")
include("functions.jl")
include("twopointfluxfvmsystem.jl")
export FVMGraph
export TwoPointFluxFVMSystem
export FVMParameters
export FVMNewtonControl
export unknowns
export fbernoulli
export fbernoulli_pm
export solve
export Dirichlet

end

