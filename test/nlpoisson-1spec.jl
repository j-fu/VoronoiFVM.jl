using Printf
using TwoPointFluxFVM

if !isinteractive()
    using PyPlot
end


mutable struct MyParameters1 <:FVMParameters
    number_of_species::Int64
    eps::Float64 
    param::Float64
    function MyParameters1()
        new(1,1,1)
    end
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
    
    
    function reaction!(this::MyParameters,f,u)
        f[1]=u[1]^2
    end
    
    function flux!(this::MyParameters,f,uk,ul)
        f[1]=this.eps*(uk[1]^2-ul[1]^2)
    end 
    
    function source!(this::MyParameters,f,x)
        f[1]=1.0e-4*x
    end 
    
    
    parameters.param=1.0e-5
    
    sys=TwoPointFluxFVMSystem(geom,parameters=parameters, flux=flux!, reaction=reaction!)
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
            plot(geom.Nodes[1,:],U[1,:])
            pause(1.0e-10)
        end
    end
end


if !isinteractive()
    @time run_1spec(n=100,pyplot=true)
    waitforbuttonpress()
end



