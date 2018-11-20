using Printf
using TwoPointFluxFVM
if isinteractive()
    using PyPlot
end


mutable struct Test2DPhysics <:FVMPhysics
    @AddFVMPhysicsBaseClassFields
    eps::Float64 
    Test2DPhysics()=Test2DPhysics(new())
end

function reaction!(this::Test2DPhysics,f,u)
    f[1]=u[1]^2
end

function flux!(this::Test2DPhysics,f,uk,ul)
    f[1]=this.eps*(uk[1]^2-ul[1]^2)
end 

function source!(this::Test2DPhysics,f,x)
    x1=x[1]-0.5
    x2=x[2]-0.5
    f[1]=exp(-20*(x1^2+x2^2))
end 
    

function Test2DPhysics(this)
    FVMPhysicsBase(this,1)
    this.eps=1
    this.reaction=reaction!
    this.flux=flux!
    this.source=source!
    return this
end


function run_test2d(;n=10,pyplot=false,verbose=false)
    
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)


    geom=FVMGraph(X,Y)
    
    physics=Test2DPhysics()
    
   
    

    
    sys=TwoPointFluxFVMSystem(geom,physics)
    sys.boundary_values[1,2]=0.1
    sys.boundary_values[1,4]=0.1
    
    sys.boundary_factors[1,2]=Dirichlet
    sys.boundary_factors[1,4]=Dirichlet
    
    inival=unknowns(sys)
    inival.=0.5

    physics.eps=1.0e-2

    control=FVMNewtonControl()
    control.verbose=verbose
    control.tol_linear=1.0e-5
    control.max_lureuse=10
    tstep=0.01
    time=0.0
    u15=0
    while time<1.0
        time=time+tstep
        U=solve(sys,inival,control=control,tstep=tstep)
        u15=U[15]
        for i in eachindex(U)
            inival[i]=U[i]
        end
        if verbose
            @printf("time=%g\n",time)
        end

        tstep*=1.0
        if pyplot
            levels=collect(0:0.01:1)
            PyPlot.clf()
            contourf(X,Y,reshape(U,length(X),length(Y)), cmap=ColorMap("hot"),levels=levels)
            colorbar()
            pause(1.0e-10)
        end
    end
    return u15
end

