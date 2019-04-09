module TwoSpeciesNonlinearPoisson

using Printf
using VoronoiFVM



if isinteractive()
    using PyPlot
end

mutable struct Physics   <: VoronoiFVM.Physics
    flux::Function
    source::Function
    reaction::Function
    storage::Function
    eps::Array{Float64,1}
    Physics()=new()
end


function main(;n=100,pyplot=false,verbose=false,dense=false)
    h=1/n
    grid=VoronoiFVM.Grid(collect(0:h:1))
    
        
    physics=Physics()
    physics.eps=[1,1]


    physics.reaction=function(physics,node,f,u)
        f[1]=u[1]*u[2]
        f[2]=-u[1]*u[2]
    end

    physics.flux=function(physics,edge,f,uk,ul)   
        f[1]=physics.eps[1]*(uk[1]-ul[1])*(0.01+uk[2]+ul[2])
        f[2]=physics.eps[2]*(uk[2]-ul[2])*(0.01+uk[1]+ul[1])
    end 
    
    physics.source=function(physics,node,f)
        f[1]=1.0e-4*(0.01+node.coord[1])
        f[2]=1.0e-4*(0.01+1.0-node.coord[1])
    end
    
    physics.storage=function(physics,node, f,u)
        f[1]=u[1]
        f[2]=u[2]
    end
    

    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics,2)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics,2)
    end
    
    add_species(sys,1,[1])
    add_species(sys,2,[1])
    
    sys.boundary_values[1,1]=1.0
    sys.boundary_values[1,2]=0.0
    
    sys.boundary_factors[1,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[1,2]=VoronoiFVM.Dirichlet
    
    sys.boundary_values[2,1]=1.0
    sys.boundary_values[2,2]=0.0
    
    sys.boundary_factors[2,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[2,2]=VoronoiFVM.Dirichlet
    
    inival=unknowns(sys)
    inival.=0
    
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    control.damp_initial=0.1
    u5=0
    for eps in [1.0,0.1,0.01]
        physics.eps=[eps,eps]
        U=solve(sys,inival,control=control)
        inival.=U
        if pyplot
            clf()
            plot(grid.coord[1,:],U[1,:])
            plot(grid.coord[1,:],U[2,:])
            pause(1.0e-10)
        end
        u5=U[5]
    end
    return u5
end


    
end
