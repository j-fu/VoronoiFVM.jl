#=

# Comparison with DifferentialEquations.jl
([source code](SOURCE_URL))

Solve the nonlinear diffusion equation

```math
\partial_t u -\Delta u^m = 0
```
in $\Omega=(-1,1)$ with homogeneous Neumann boundary conditons using the implicit Euler method.

This equation is also called  "porous medium equation". 
The Barenblatt solution is an exact solution of this problem which for m>1 has a finite support.
We initialize this problem with the exact solution for $t=t_0=0.001$.

(see Barenblatt, G. I. "On nonsteady motions of gas and fluid in porous medium." Appl. Math. and Mech.(PMM) 16.1 (1952): 67-78.)

Here, we compare the implicit Euler approach in VoronoiFVM with the ODE solvers in [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) and demonstrate the possibility to use VoronoiFVM to define differential
operators compatible with its [ODEFunction](https://diffeq.sciml.ai/stable/features/performance_overloads/#ODEFunction)
interface.

Note that this example requires  PyPlot and DifferentialEquations to be installed.
=#


module TestDiffEq

using VoronoiFVM
using DifferentialEquations
using LinearAlgebra
using Printf
using PyPlot

import Base:push!





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


function create_porous_medium_problem(n,m,unknown_storage)
    h=1.0/convert(Float64,n/2)
    X=collect(-1:h:1)
    grid=VoronoiFVM.Grid(X)

    function flux!(f,u0,edge)
        u=unknowns(edge,u0)
        f[1]=u[1,1]^m-u[1,2]^m
    end

    storage!(f,u,node)= f[1]=u[1]
    
    physics=VoronoiFVM.Physics(flux=flux!,storage=storage!)
    
    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)
    enable_species!(sys,1,[1])
    sys,X
end


function run_vfvm(;n=20,m=2,t0=0.001, tend=0.01,tstep=1.0e-6,unknown_storage=:dense)
    
    sys,X=create_porous_medium_problem(n,m,unknown_storage)

    inival=unknowns(sys)
    inival[1,:].=map(x->barenblatt(x,t0,m),X)

    solution=unknowns(sys)
    control=VoronoiFVM.NewtonControl()
    control.verbose=false
    control.Δt=tstep
    control.Δu_opt=0.05
    control.Δt_min=tstep

    times=collect(t0:t0:tend)
    times=[t0,tend]
    sol=VoronoiFVM.solve(inival,sys,times,control=control,store_all=true)
    err=norm(sol[1,:,end]-map(x->barenblatt(x,tend,m),X))
    sol,X,err
end



function run_diffeq(;n=20,m=2, t0=0.001,tend=0.01, unknown_storage=:dense,solver=nothing)
    sys,X=create_porous_medium_problem(n,m,unknown_storage)
    inival=unknowns(sys)
    inival[1,:].=map(x->barenblatt(x,t0,m),X)
    tspan = (t0,tend)
    sol=VoronoiFVM.solve(DifferentialEquations,inival,sys,tspan,solver=solver)
    err=norm(sol[1,:,end]-map(x->barenblatt(x,tend,m),X))
    sol, X,err
end


function main(;m=2,n=20, solver=nothing, unknown_storage=:dense)
    function plotsol(sol,X)
        f=sol[1,:,:]'
        contourf(X,sol.t,f,0:0.1:10,cmap=:summer)
        contour(X,sol.t,f,0:1:10,colors=:black)
    end

    clf()
    subplot(121)
    
    t=@elapsed begin
        sol,X,err=run_vfvm(m=m,n=n, unknown_storage=unknown_storage)
    end
    title(@sprintf("VoronoiFVM: %.0f ms e=%.2e",t*1000,err))
    plotsol(sol,X)
    
    subplot(122)
    t=@elapsed begin
        sol,X,err=run_diffeq(m=m,n=n,solver=solver, unknown_storage=unknown_storage)
    end
    plotsol(sol,X)
    title(@sprintf("DifferentialEq: %.0f ms, e=%.2e",t*1000,err))
    
    gcf().set_size_inches(8,4)
    gcf()
end
end
