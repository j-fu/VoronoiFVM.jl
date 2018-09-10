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




mutable struct MyParameters <:FVMParameters
    number_of_species::Int64
    eps::Array{Float64,1}
    MyParameters()=new(2)
end


function reaction!(this::MyParameters,f,u)
    f[1]=u[1]*u[2]
    f[2]=-u[1]*u[2]
end

function flux!(this::MyParameters,f,uk,ul)   
    f[1]=this.eps[1]*(uk[1]-ul[1])*(0.01+uk[2]+ul[2])
    f[2]=this.eps[2]*(uk[2]-ul[2])*(0.01+uk[1]+ul[1])
end 

function source!(this::MyParameters,f,x)
    f[1]=1.0e-4*(0.01+x[1])
    f[2]=1.0e-4*(0.01+1.0-x[1])
end 



parameters=MyParameters()


sys=TwoPointFluxFVMSystem(geom,parameters=parameters,
                          flux=flux!,
                          reaction=reaction!,
                          source=source!)
                        
sys.dirichlet_values[1,1]=1.0
sys.dirichlet_values[1,2]=0.0
                       
sys.dirichlet_values[2,1]=0.0
sys.dirichlet_values[2,2]=1.0

inival=unknowns(sys)
inival.=0


for eps in [1.0,0.1,0.01]
    parameters.eps=[eps,eps]
    U=solve(sys,inival)
    if args["pyplot"]
        plot(geom.Points[1,:],U[1,:])
        plot(geom.Points[1,:],U[2,:])
        pause(1.0e-10)
        waitforbuttonpress()
    end
end

