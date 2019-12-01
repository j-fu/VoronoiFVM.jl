#=

# 106: 1D Nonlinear Diffusion equation

Solve the nonlinear diffusion equation

```math
\partial_t u -\Delta u^m = 0
```
in $\Omega=(-1,1)$ with boundary condition $u(-1)=0$ and $u(1)=0$ using the implicit Euler method.

This equation is also called  "porous medium equation". 
The Barenblatt solution is an exact solution of this problem which for m>1 has a finite support.
We initialize this problem with the exact solution for $t=t_0=0.001$.

(see Barenblatt, G. I. "On nonsteady motions of gas and fluid in porous medium." Appl. Math. and Mech.(PMM) 16.1 (1952): 67-78.)
=# 

module Example106_NonlinearDiffusion1D
using Printf
using VoronoiFVM


function barenblatt(x,t,m)
    tx=t^(-1.0/(m+1.0))
    xx=x*tx
    xx=xx*xx
    xx=1- xx*(m-1)/(2.0*m*(m+1));
    if xx<0.0
        xx=0.0
    end
    return tx*xx^(1.0/(m-1.0))
end


function main(;n=20,m=2,Plotter=nothing,verbose=false, dense=false,tend=0.01,tstep=0.0001)
    
    ## Create a one-dimensional discretization
    h=1.0/convert(Float64,n/2)
    X=collect(-1:h:1)
    grid=VoronoiFVM.Grid(X)

    ## Flux function which describes the flux
    ## between neigboring control volumes
    function flux!(f,u,edge,data)
        uk=viewK(edge,u)  
        ul=viewL(edge,u)
        f[1]=uk[1]^m-ul[1]^m
    end

    ## Storage term
    function storage!(f,u,node,data)
        f[1]=u[1]
    end
    
    ## Create a physics structure
    physics=VoronoiFVM.Physics(
        flux=flux!,
        storage=storage!)
    
    
    ## Create a finite volume system - either
    ## in the dense or  the sparse version.
    ## The difference is in the way the solution object
    ## is stored - as dense or as sparse matrix
    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics)
    end

    ## Add species 1 to region 1
    enable_species!(sys,1,[1])
    
    ## Set boundary conditions
    sys.boundary_values[1,1]=0.0
    sys.boundary_values[1,2]=0.0
    sys.boundary_factors[1,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[1,2]=VoronoiFVM.Dirichlet

    ## Create a solution array
    inival=unknowns(sys)
    solution=unknowns(sys)
    t0=0.001
    ## Broadcast the initial value
    inival[1,:].=map(x->barenblatt(x,t0,m),X)


    ## Create solver control info
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    time=t0
    while time<tend
        time=time+tstep
        solve!(solution,inival,sys,control=control,tstep=tstep)
        inival.=solution
        if verbose
            @printf("time=%g\n",time)
        end
        if isplots(Plotter)
            p=Plotter.plot(X,
                         solution[1,:],
                         label="numerical",
                         title=@sprintf("Nonlinear Diffusion t=%.5f",time),
                         grid=true)
            Plotter.plot!(p,X,
                        map(x->barenblatt(x,time,m),X),
                        label="exact",
                        show=true)
        end
    end
    return sum(solution)
end


function test()
    testval=46.66666666647518
    main(dense=false) ≈ testval && main(dense=true) ≈ testval
end

# End of module
end 

