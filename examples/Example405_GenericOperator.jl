#=
# 405: Generic Operator: 1D Nonlinear Poisson equation 
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

module Example405_GenericOperator
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize

function main(;n=10,Plotter=nothing,verbose=false, unknown_storage=:sparse)
    
    ## Create a one-dimensional discretization
    h=1.0/convert(Float64,n)
    X=collect(0:h:1)
    grid=VoronoiFVM.Grid(X)

    ## A parameter which is "passed" to the flux function via scope
    ϵ=1.0e-3

    ## This generic operator works on the full solution seen as linear vector, and indexing
    ## into it needs to be performed with the help of idx (defined below for a solution vector)
    ## Here, instead of the flux function we provide a "generic operator"
    ## which provides the stiffness part of the problem. Its sparsity is detected automatically
    ## using SparsityDetection.jl 
    function generic_operator!(f,u,sys)
        f.=0.0
        for i=1:length(X)-1
            du=ϵ*(u[idx[1,i]]-u[idx[1,i+1]])/(X[i+1]-X[i])
            f[idx[1,i]]+=du
            f[idx[1,i+1]]-=du
        end
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
        generic=generic_operator!,
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
    
    idx=unknown_indices(solution)
    ## Create solver control info
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose

    ## Stationary solution of the problem
    solve!(solution,inival,sys, control=control)
    
    scalarplot(grid,solution[1,:],title="Nonlinear Poisson",Plotter=Plotter)
    return sum(solution)
end

function test()
    testval=1.5247901344230088
    main(unknown_storage=:sparse) ≈ testval && main(unknown_storage=:dense) ≈ testval
end

end 

