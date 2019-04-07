module TwoSpeciesTestFunctions

using Printf
using TwoPointFluxFVM



if isinteractive()
    using PyPlot
end

mutable struct Physics   <: TwoPointFluxFVM.Physics
    flux::Function
    source::Function
    reaction::Function
    storage::Function
    eps::Array{Float64,1}
    Physics()=new()
end


function main(;n=100,pyplot=false,verbose=false,dense=false)
    h=1/n
    grid=TwoPointFluxFVM.Grid(collect(0:h:1))
    
        
    physics=Physics()
    physics.eps=[1,1.0e-1]


    physics.reaction=function(physics,node,f,u)
        f[1]=10*(u[1]-u[2])
        f[2]=10*(u[2]-u[1])
    end

    physics.flux=function(physics,edge,f,uk,ul)   
        f[1]=physics.eps[1]*(uk[1]-ul[1])
        f[2]=physics.eps[2]*(uk[2]-ul[2])
    end 
    
    
    physics.storage=function(physics,node, f,u)
        f[1]=u[1]
        f[2]=u[2]
    end
    

    if dense
        sys=TwoPointFluxFVM.DenseSystem(grid,physics,2)
    else
        sys=TwoPointFluxFVM.SparseSystem(grid,physics,2)
    end
    
    add_species(sys,1,[1])
    add_species(sys,2,[1])
    
    sys.boundary_values[1,1]=0.01
    sys.boundary_values[1,2]=0.0
    
    sys.boundary_factors[1,1]=0
    sys.boundary_factors[1,2]=0
    
    sys.boundary_values[2,1]=0.0
    sys.boundary_values[2,2]=0.0
    
    sys.boundary_factors[2,1]=0
    sys.boundary_factors[2,2]=TwoPointFluxFVM.Dirichlet

    factory=TwoPointFluxFVM.TestFunctionFactory(sys)
    tf1=testfunction(factory,[2],[1])
    tf2=testfunction(factory,[1],[2])
    
    if pyplot
        clf()
        plot(grid.coord[1,:],tf1)
        pause(1.0e-10)
        waitforbuttonpress()
        plot(grid.coord[1,:],tf2)
        pause(1.0e-10)
        waitforbuttonpress()
    end
    
    inival=unknowns(sys)
    inival[2,:].=0.1
    inival[1,:].=0.1
    
    control=TwoPointFluxFVM.NewtonControl()
    control.verbose=verbose
    control.damp_initial=0.1
    u5=0
    for eps in [1.0,0.1,0.01]
        physics.eps=[eps,eps]
        U=solve(sys,inival,control=control)
        println(integrate(sys,tf1,U))
        println(integrate(sys,tf2,U))
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
