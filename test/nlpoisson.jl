# Command line argument parser, get it via Pkg.add
using ArgParse
using TwoPointFluxFVM

argdef = ArgParseSettings()
add_arg_table(argdef,
    "--pyplot", Dict(:help => "call python visualization",
#                    :default => false,
                    :action => :store_true),
    "--n", Dict(:help => "problem size",
                :arg_type => Int,
                :default => 11)

)

args = parse_args(argdef)

if args["pyplot"]
    using PyPlot
end



geom=FVMGraph(collect(0:0.01:1))

#
# Structure containing "physics" information
#
#  - \nabla \eps  \nabla u + reaction(u) = source(x)
#




mutable struct MyPhysics <:FVMPhysics
  eps::Float64 
  param::Float64
  nspec::Function
  reaction::Function
  flux::Function
  source::Function
  MyPhysics()=new()
end


physics=MyPhysics()
physics.nspec= function(this) return 1 end


physics.reaction=function(this::MyPhysics,f,u)
   f[1]=u[1]^2+this.param*u[1]
end

physics.flux=function(this::MyPhysics,f,uk,ul)
    f[1]=this.eps*(uk[1]^2-ul[1]^2)
end 

physics.source=function(this::MyPhysics,f,x)
    f[1]=1.0e-4*x[1]
end 




physics.param=1.0e-5

sys=TwoPointFluxFVMSystem(geom,physics)
sys.dirichlet_values[1,1]=1.0
sys.dirichlet_values[2,1]=0.0

inival=unknowns(sys)
inival.=0

for eps in [1.0,0.1,0.01,0.0001]
    physics.eps=eps
    U=solve(sys,inival)
    if args["pyplot"]
        plot(geom.Points[:,1],U[:,1])
        pause(1.0e-10)
        waitforbuttonpress()
    end
end
