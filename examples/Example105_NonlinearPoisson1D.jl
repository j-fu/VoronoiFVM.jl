#=
# 105: 1D Nonlinear Poisson equation 
([source code](SOURCE_URL))

Solve the nonlinear Poisson equation

```math
-\nabla \varepsilon \nabla u + e^{u}-e^{-u} = f
```
in $\Omega=(0,1)$ with boundary condition $u(0)=0$ and $u(1)=1$ with 
```math
f(x)=
    \begin{cases}
    1&,x>0.5\\
    -1&, x<0.5
    \end{cases}.
```
    
This stationary problem is an example of a nonlinear Poisson equation or Poisson-Boltzmann equation.
Such equation occur e.g. in simulations of electrochemical systems and semicondutor devices.
 
=#

module Example105_NonlinearPoisson1D
using Printf
using VoronoiFVM
using ExtendableGrids
using .GridVisualize

function main(;n=10,Plotter=nothing,verbose=false, unknown_storage=:sparse)
    
    ## Create a one-dimensional discretization
    h=1.0/convert(Float64,n)
    grid=VoronoiFVM.Grid(collect(0:h:1))

    ## A parameter which is "passed" to the flux function via scope
    ϵ=1.0e-3

    ## Flux function which describes the flux
    ## between neigboring control volumes
    function flux!(f,u0,edge)
        u=unknowns(edge,u0)
        f[1]=ϵ*(u[1,1]-u[1,2])
    end

    ## Source term
    function source!(f,node)
        if node[1]<=0.5
            f[1]=1
        else
            f[1]=-1
        end
    end
    
    ## Reaction term
    function reaction!(f,u,node)
        f[1]=exp(u[1]) - exp(-u[1]) 
    end
    
    ## Create a physics structure
    physics=VoronoiFVM.Physics(
        flux=flux!,
        source=source!,
        reaction=reaction!)
    

    ## Create a finite volume system - either
    ## in the dense or  the sparse version.
    ## The difference is in the way the solution object
    ## is stored - as dense or as sparse matrix

    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)

    ## Add species 1 to region 1
    enable_species!(sys,1,[1])

    ## Set boundary conditions
    boundary_dirichlet!(sys,1,1,0.0)
    boundary_dirichlet!(sys,1,2,1.0)

    ## Create a solution array
    inival=unknowns(sys,inival=0.5)
    solution=unknowns(sys)

    ## Create solver control info
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose

    ## Stationary solution of the problem
    solve!(solution,inival,sys, control=control)

    visualize(grid,solution[1,:],title="Nonlinear Poisson", Plotter=Plotter)

    return sum(solution)
end

function test()
    testval=1.5247901344230088
    main(unknown_storage=:sparse) ≈ testval && main(unknown_storage=:dense) ≈ testval
end

end 

