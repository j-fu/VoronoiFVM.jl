module TwoSpeciesNonlinearPoisson

using Printf
using TwoPointFluxFVM

if isinteractive()
    using PyPlot
end

mutable struct Physics <:FVMPhysics
    @AddFVMPhysicsBaseClassFields
    eps::Array{Float64,1}
    Physics()=Physics(new())
end

function reaction!(this::Physics,f,u)
    f[1]=u[1]*u[2]
    f[2]=-u[1]*u[2]
end

function flux!(this::Physics,f,uk,ul)   
    f[1]=this.eps[1]*(uk[1]-ul[1])*(0.01+uk[2]+ul[2])
    f[2]=this.eps[2]*(uk[2]-ul[2])*(0.01+uk[1]+ul[1])
end 

function source!(this::Physics,f,x)
    f[1]=1.0e-4*(0.01+x[1])
    f[2]=1.0e-4*(0.01+1.0-x[1])
end 
    

function Physics(this::Physics)
    FVMPhysicsBase(this,2)
    this.eps=[1,1]
    this.flux=flux!
    this.reaction=reaction!
    this.source=source!
    return this
end

function main(;n=10,pyplot=false,verbose=true)
    
    geom=FVMGraph(collect(0:0.01:1))
    
        
    physics=Physics()
    
    
    sys=TwoPointFluxFVMSystem(geom,physics)
    
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
    
    control=FVMNewtonControl()

    u5=0
    for eps in [1.0,0.1,0.01]
        physics.eps=[eps,eps]
        U0=solve(sys,inival)
        U=bulk_unknowns(sys,U0)
        if pyplot
            clf()
            plot(geom.node_coordinates[1,:],U[1,:])
            plot(geom.node_coordinates[1,:],U[2,:])
            pause(1.0e-10)
        end
        u5=U0[5]
    end
    return u5
end


    
end
