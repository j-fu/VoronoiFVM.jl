using Printf
using TwoPointFluxFVM

if isinteractive()
    using PyPlot
end

mutable struct NLPoisson1SpecPhysics <:FVMPhysics
    @AddDefaultFVMPhysics
    eps::Float64 
    param::Float64
    NLPoisson1SpecPhysics()=NLPoisson1SpecPhysics(new())
end


function reaction!(this::NLPoisson1SpecPhysics,f,u)
    f[1]=u[1]^2
end

function flux!(this::NLPoisson1SpecPhysics,f,uk,ul)
    f[1]=this.eps*(uk[1]^2-ul[1]^2)
end 

function source!(this::NLPoisson1SpecPhysics,f,x)
    f[1]=1.0e-4*x[1]
end 


function NLPoisson1SpecPhysics(this::NLPoisson1SpecPhysics)
    DefaultFVMPhysics(this,1)
    this.eps=1
    this.param=1
    this.flux=flux!
    this.reaction=reaction!
    this.source=source!
    return this
end


function run_nlpoisson_1spec(;n=10,pyplot=false,verbose=false)
    h=1.0/convert(Float64,n)
    geom=FVMGraph(collect(0:h:1))
    
    #
    # Structure containing "physics" information
    #
    #  - \nabla \eps  \nabla u + reaction(u) = source(x)
    #
    
    
    physics=NLPoisson1SpecPhysics()
    physics.param=1.0e-5

    sys=TwoPointFluxFVMSystem(geom,physics)
    sys.boundary_values[1,1]=1.0
    sys.boundary_values[1,2]=0.5
    
    sys.boundary_factors[1,1]=Dirichlet
    sys.boundary_factors[1,2]=Dirichlet
    
    inival=unknowns(sys)
    inival.=0.5

    physics.eps=1.0e-2

    control=FVMNewtonControl()
    control.verbose=verbose
    tstep=1.0e-2
    times=collect(0.0:tstep:1.0)
    u5=0
    for it=2:length(times)
        U=solve(sys,inival,control=control,tstep=tstep)
        u5=U[5]
        for i in eachindex(U)
            inival[i]=U[i]
        end
        if verbose
            @printf("time=%g\n",times[it])
        end
        if pyplot
            PyPlot.clf()
            plot(geom.node_coordinates[1,:],U)
            pause(1.0e-10)
        end
    end
    return u5
end
