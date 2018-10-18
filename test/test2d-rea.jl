using Printf
using TwoPointFluxFVM
using PyPlot


mutable struct Test2DReaPhysics <:FVMPhysics
    @AddDefaultFVMPhysics
    k::Float64
    eps::Float64 
    Test2DReaPhysics()=Test2DReaPhysics(new())
end


function reaction!(this::Test2DReaPhysics,f,u)
    f[1]=this.k*(u[1]-u[2])
    f[2]=this.k*(u[2]-u[1])
end

function flux!(this::Test2DReaPhysics,f,uk,ul)
    f[1]=this.eps*(uk[1]-ul[1])
    f[2]=this.eps*(uk[2]-ul[2])
end 

function source!(this::Test2DReaPhysics,f,x)
    x1=x[1]-0.5
    x2=x[2]-0.5
    f[1]=exp(-20*(x1^2+x2^2))
end 



function Test2DReaPhysics(this)
    DefaultFVMPhysics(this,2)
    this.eps=1
    this.k=1.0
    this.reaction=reaction!
    this.flux=flux!
    this.source=source!
    return this
end


function run_test2d_rea(;n=10,pyplot=false,verbose=true)
    
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)
    


    geom=FVMGraph(X,Y)

    
    physics=Test2DReaPhysics()
    
   
    

    
    sys=TwoPointFluxFVMSystem(geom,physics)
    # sys.boundary_values[2,2]=0
    # sys.boundary_values[2,4]=0
    
    # sys.boundary_factors[2,2]=Dirichlet
    # sys.boundary_factors[2,4]=Dirichlet
    
    inival=unknowns(sys)
    inival.=0.0
    
    physics.eps=1.0e-2
    
    control=FVMNewtonControl()
    control.verbose=verbose
    control.tol_linear=1.0e-5
    control.max_lureuse=0
    tstep=0.01
    time=0.0
    istep=0
    u15=0
    while time<1
        time=time+tstep
        U=solve(sys,inival,control=control,tstep=tstep)
        inival.=U
        # for i in eachindex(U)
        #     inival[i]=U[i]
        # end
        if verbose
            @printf("time=%g\n",time)
        end
        u15=U[15]
        tstep*=1.0
        istep=istep+1
        if pyplot && istep%10 == 0
            U_bulk=bulk_unknowns(sys,U)
            if verbose
                @printf("max1=%g max2=%g\n",maximum(U_bulk[1,:]),maximum(U_bulk[2,:]))
            end
            levels1=collect(0:0.01:1)
            PyPlot.clf()
            subplot(211)
            contourf(X,Y,reshape(U_bulk[1,:],length(X),length(Y)), cmap=ColorMap("hot"))
            colorbar()
            levels2=collect(0:0.001:0.1)
            subplot(212)
            contourf(X,Y,reshape(U_bulk[2,:],length(X),length(Y)), cmap=ColorMap("hot"))
            colorbar()

            pause(1.0e-10)
        end
    end
    return u15
end

