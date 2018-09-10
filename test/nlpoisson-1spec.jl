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



geom=FVMGraph(collect(0:0.1:1))

#
# Structure containing "physics" information
#
#  - \nabla \eps  \nabla u + reaction(u) = source(x)
#




mutable struct MyParameters <:FVMParameters
    number_of_species::Int64
    eps::Float64 
    param::Float64
    function MyParameters()
        new(1,1,1)
    end
end


parameters=MyParameters()


function reaction!(this::MyParameters,f,u)
    f[1]=u[1]^2+this.param*u[1]
end

function flux!(this::MyParameters,f,uk,ul)
    f[1]=this.eps*(uk[1]^2-ul[1]^2)
end 

function source!(this::MyParameters,f,x)
    f[1]=1.0e-4*x[1]
end 


parameters.param=1.0e-5

sys=TwoPointFluxFVMSystem(geom,parameters=parameters, flux=flux!, reaction=reaction!)
sys.dirichlet_values[1,1]=1.0
sys.dirichlet_values[1,2]=0.0

inival=unknowns(sys)
inival.=1

for eps in [1.0,0.1,0.01]
    parameters.eps=eps
    U=solve(sys,inival)
    if args["pyplot"]
        plot(geom.Points[1,:],U[1,:])
        pause(1.0e-10)
        waitforbuttonpress()
    end
end
