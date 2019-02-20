module TwoSpeciesNonlinearPoisson

using Printf
using TwoPointFluxFVM
const Node=TwoPointFluxFVM.Node
const Edge=TwoPointFluxFVM.Edge

if isinteractive()
    using PyPlot
end

mutable struct Physics <:TwoPointFluxFVM.Physics
    TwoPointFluxFVM.@AddPhysicsBaseClassFields
    eps::Array{Float64,1}
    Physics()=Physics(new())
end

function reaction!(this::Physics,node::Node,f,u)
    f[1]=u[1]*u[2]
    f[2]=-u[1]*u[2]
end

function flux!(this::Physics,edge::Edge,f,uk,ul)   
    f[1]=this.eps[1]*(uk[1]-ul[1])*(0.01+uk[2]+ul[2])
    f[2]=this.eps[2]*(uk[2]-ul[2])*(0.01+uk[1]+ul[1])
end 

function source!(this::Physics,node::Node,f)
    f[1]=1.0e-4*(0.01+node.coord[1])
    f[2]=1.0e-4*(0.01+1.0-node.coord[1])
end 
    

function Physics(this::Physics)
    TwoPointFluxFVM.PhysicsBase(this,2)
    this.eps=[1,1]
    this.flux=flux!
    this.reaction=reaction!
    this.source=source!
    return this
end

function main(;n=10,pyplot=false,verbose=true)
    
    geom=TwoPointFluxFVM.Graph(collect(0:0.01:1))
    
        
    physics=Physics()
    
    
    sys=TwoPointFluxFVM.System(geom,physics)
    
    sys.boundary_values[1,1]=1.0
    sys.boundary_values[1,2]=0.0
    
    sys.boundary_factors[1,1]=TwoPointFluxFVM.Dirichlet
    sys.boundary_factors[1,2]=TwoPointFluxFVM.Dirichlet
    
    sys.boundary_values[2,1]=1.0
    sys.boundary_values[2,2]=0.0
    
    sys.boundary_factors[2,1]=TwoPointFluxFVM.Dirichlet
    sys.boundary_factors[2,2]=TwoPointFluxFVM.Dirichlet
    
    inival=unknowns(sys)
    inival.=0
    
    control=TwoPointFluxFVM.NewtonControl()

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
