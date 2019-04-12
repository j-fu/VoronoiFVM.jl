module TwoSpeciesNonlinearPoisson

using Printf
using VoronoiFVM



if isinteractive()
    using PyPlot
end



function main(;n=100,pyplot=false,verbose=false,dense=false)
    h=1/n
    grid=VoronoiFVM.Grid(collect(0:h:1))
    
    
    eps=[1.0,1.0]
    
    physics=VoronoiFVM.Physics(num_species=2,
                               
                               reaction=function(f,u,node,data)
                               f[1]=u[1]*u[2]
                               f[2]=-u[1]*u[2]
                               end,
                               
                               flux=function(f,u,edge,data)   
                               nspecies=2
                               uk=viewK(2,u)
                               ul=viewL(2,u)
                               f[1]=eps[1]*(uk[1]-ul[1])*(0.01+uk[2]+ul[2])
                               f[2]=eps[2]*(uk[2]-ul[2])*(0.01+uk[1]+ul[1])
                               end,
                               
                               source=function(f,node,data)
                               f[1]=1.0e-4*(0.01+node.coord[1])
                               f[2]=1.0e-4*(0.01+1.0-node.coord[1])
                               end,
                               
                               storage=function(f,u,node,data)
                               f[1]=u[1]
                               f[2]=u[2]
                               end
                               )
    
    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics)
    end
    
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])
    
    sys.boundary_values[1,1]=1.0
    sys.boundary_values[1,2]=0.0
    
    sys.boundary_factors[1,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[1,2]=VoronoiFVM.Dirichlet
    
    sys.boundary_values[2,1]=1.0
    sys.boundary_values[2,2]=0.0
    
    sys.boundary_factors[2,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[2,2]=VoronoiFVM.Dirichlet
    
    inival=unknowns(sys)
    U=unknowns(sys)
    inival.=0
    
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    control.damp_initial=0.1
    u5=0
    for xeps in [1.0,0.1,0.01]
        eps=[xeps,xeps]
        solve!(U,inival,sys,control=control)
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
