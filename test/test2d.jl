using Printf
using TwoPointFluxFVM

if !isinteractive()
    using PyPlot
end


mutable struct Test2DParameters <:FVMParameters
    number_of_species::Int64
    eps::Float64 
    function Test2DParameters()
        new(1,1)
    end
end

function run_test2d(;n=100,pyplot=false)
    
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)


    geom=FVMGraph(X,Y)
    
    parameters=Test2DParameters()
    
   
    
    function reaction!(this::Test2DParameters,f,u)
        f[1]=u[1]^2
    end
    
    function flux!(this::Test2DParameters,f,uk,ul)
        f[1]=this.eps*(uk[1]^2-ul[1]^2)
    end 
    
    function source!(this::Test2DParameters,f,x)
        x1=x[1]-0.5
        x2=x[2]-0.5
        f[1]=exp(-20*(x1^2+x2^2))
    end 
    
    
    sys=TwoPointFluxFVMSystem(geom,parameters=parameters, flux=flux!, reaction=reaction!,source=source!)
    sys.boundary_values[1,2]=0.1
    sys.boundary_values[1,4]=0.1
    
    sys.boundary_factors[1,2]=Dirichlet
    sys.boundary_factors[1,4]=Dirichlet
    
    inival=unknowns(sys)
    inival.=0.5

    parameters.eps=1.0e-2

    control=FVMNewtonControl()
    control.verbose=true
    control.lin_tolerance=1.0e-5
    control.max_lureuse=10
    tstep=0.01
    time=0.0
    while time<1.0
        time=time+tstep
        U=solve(sys,inival,control=control,tstep=tstep)
        for i in eachindex(U)
            inival[i]=U[i]
        end
        @printf("time=%g\n",time)

        tstep*=1.0
        @time if pyplot
            levels=collect(0:0.01:1)
            PyPlot.clf()
            contourf(X,Y,reshape(U[1,:],length(X),length(Y)), cmap=ColorMap("hot"),levels=levels)
            colorbar()
            pause(1.0e-10)
        end
    end
end


if !isinteractive()
    @time run_test2d(n=100,pyplot=true)
    waitforbuttonpress()
end



