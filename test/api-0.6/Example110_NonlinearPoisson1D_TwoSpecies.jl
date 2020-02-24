# # 110: 1D Nonlinear Poisson equation with two species
# ([source code](SOURCE_URL))
# Solve the nonlinear Poisson equation
# 
# ```math
# -\nabla (0.01+2u_2)\nabla u_1 + u_1u_2= 0.0001(0.01+x)
# ```
# 
# ```math
# -\nabla (0.01+2u_1)\nabla u_2 -+ u_1u_2 = 0.0001(1.01-x)
# ```
# 
# 
# in $\Omega=(0,1)$ with boundary condition $u_1(0)=1$, $u_2(0)=0$ and $u_1(1)=1$, $u_2(1)=1$.
# 


module Example110_NonlinearPoisson1D_TwoSpecies

using Printf
using VoronoiFVM




function main(;n=100,Plotter=nothing,verbose=false,dense=false)
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

    boundary_dirichlet!(sys,1,1,1.0)
    boundary_dirichlet!(sys,1,2,0.0)
    
    boundary_dirichlet!(sys,2,1,1.0)
    boundary_dirichlet!(sys,2,2,0.0)
    
    inival=unknowns(sys)
    U=unknowns(sys)
    inival.=0
    
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    control.damp_initial=0.1
    u5=0
    for xeps in [1.0,0.5,0.25,0.1,0.05,0.025,0.01]
        eps=[xeps,xeps]
        solve!(U,inival,sys,control=control)
        inival.=U
        if isplots(Plotter)
            p=Plotter.plot(grid.coord[1,:],U[1,:], grid=true)
            Plotter.plot!(p,grid.coord[1,:],U[2,:],show=true, title=@sprintf("\$\\varepsilon=%8.3f\$",xeps)),
            Plotter.sleep(0.2)
        end
        u5=U[5]
    end
    return u5
end

function test()
    main(dense=false) ≈ 0.7117546972922056 &&
        main(dense=true) ≈ 0.7117546972922056
end
end

