module TwoPointFluxFVM

include("fvmgraph.jl")
include("functions.jl")
include("fvmphysics.jl")
include("fvmnewtoncontrol.jl")
include("twopointfluxfvmsystem.jl")



export FVMGraph
export TwoPointFluxFVMSystem
export FVMPhysics
export FVMPhysicsBase
export FVMNewtonControl
export unknowns
export bulk_unknowns
export boundary_unknowns
export fbernoulli
export fbernoulli_pm
export solve
export integrate
export Dirichlet
export @AddFVMPhysicsBaseClassFields
end

