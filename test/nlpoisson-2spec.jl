using Printf
using TwoPointFluxFVM

if !isinteractive()
    using PyPlot
end



mutable struct MyParameters2 <:FVMParameters
    number_of_species::Int64
    eps::Array{Float64,1}
    MyParameters2()=new(2)
end


function run_2spec(;n=100,pyplot=false)

    geom=FVMGraph(collect(0:0.01:1))

    
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
    
    
    
    parameters=MyParameters2()
    
    
    sys=TwoPointFluxFVMSystem(geom,parameters=parameters,
                              flux=flux!,
                              reaction=reaction!,
                              source=source!)
    
    
    
    sys.boundary_values[1,1]=1.0
    sys.boundary_values[1,2]=0.0
    
    sys.boundary_factors[1,1]=Dirichlet
    sys.boundary_factors[1,2]=Dirichlet
    
    sys.boundary_values[2,1]=1.0
    sys.boundary_values[2,2]=0.0
    
    sys.boundary_factors[2,1]=Dirichlet
    sys.boundary_factors[2,2]=Dirichlet
    
    inival=unknowns(sys)
    inival.=0
    
    
    for eps in [1.0,0.1,0.01]
        parameters.eps=[eps,eps]
        U=solve(sys,inival)
        if pyplot
            plot(geom.Nodes[1,:],U[1,:])
            plot(geom.Nodes[1,:],U[2,:])
            pause(1.0e-10)
        end
    end
end

if !isinteractive()
    @time run_2spec(n=100,pyplot=true)
    waitforbuttonpress()
end

    
