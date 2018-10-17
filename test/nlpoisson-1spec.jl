using Printf
using TwoPointFluxFVM

if !isinteractive()
    using PyPlot
end


mutable struct MyParameters1 <:FVMParameters
    @AddDefaultFVMParameters
    eps::Float64 
    param::Float64
    MyParameters1()=MyParameters1(new())
end

function MyParameters1(this::MyParameters1)
    DefaultFVMParameters(this,1)
    eps=1
    param=1
    return this
end


function run_1spec(;n=100,pyplot=false)
    h=1.0/convert(Float64,n)
    geom=FVMGraph(collect(0:h:1))
    
    #
    # Structure containing "physics" information
    #
    #  - \nabla \eps  \nabla u + reaction(u) = source(x)
    #
    
    
    parameters=MyParameters1()
    
    
    function reaction!(this::MyParameters1,f,u)
        f[1]=u[1]^2
    end
    
    function flux!(this::MyParameters1,f,uk,ul)
        f[1]=this.eps*(uk[1]^2-ul[1]^2)
    end 
    
    function source!(this::MyParameters1,f,x)
        f[1]=1.0e-4*x[1]
    end 
    
    
    parameters.param=1.0e-5
    parameters.reaction=reaction!
    parameters.flux=flux!
    parameters.source=source!

    sys=TwoPointFluxFVMSystem(geom,parameters)
    sys.boundary_values[1,1]=1.0
    sys.boundary_values[1,2]=0.5
    
    sys.boundary_factors[1,1]=Dirichlet
    sys.boundary_factors[1,2]=Dirichlet
    
    inival=unknowns(sys)
    inival.=0.5

    parameters.eps=1.0e-2

    control=FVMNewtonControl()
    control.verbose=false
    tstep=1.0e-2
    times=collect(0.0:tstep:1.0)
    for it=2:length(times)
        U=solve(sys,inival,control=control,tstep=tstep)
        for i in eachindex(U)
            inival[i]=U[i]
        end
        @printf("time=%g\n",times[it])
        if pyplot
            PyPlot.clf()
            plot(geom.node_coordinates[1,:],U)
            pause(1.0e-10)
        end
    end
end


if !isinteractive()
    @time run_1spec(n=100,pyplot=true)
    waitforbuttonpress()
end



