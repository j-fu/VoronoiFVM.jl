#=

# 108: 1D Nonlinear Diffusion with OrdinaryDiffEq
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

At the moment, this code needs OrdinaryDiffEq v6.0.3.
=#


module Example108_NonlinearDiffusion1D_ODE

using VoronoiFVM
using DifferentialEquations
using LinearAlgebra
using Printf
using GridVisualize

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

    function flux!(f,u,edge)
        f[1]=u[1,1]^m-u[1,2]^m
    end

    storage!(f,u,node)= f[1]=u[1]
    sys=VoronoiFVM.System(grid,flux=flux!,storage=storage!, species=1,unknown_storage=unknown_storage)
    sys,X
end


function run_vfvm(;n=20,m=2,t0=0.001, tend=0.01,tstep=1.0e-6,unknown_storage=:dense)
    sys,X=create_porous_medium_problem(n,m,unknown_storage)
    inival=unknowns(sys)
    inival[1,:].=map(x->barenblatt(x,t0,m),X)
    sol=VoronoiFVM.solve(sys;inival,times=(t0,tend),Δt=tstep,Δu_opt=0.01,Δt_min=tstep,store_all=true,log=true)
    err=norm(sol[1,:,end]-map(x->barenblatt(x,tend,m),X))
    sol,sys,err
end


function run_diffeq(;n=20,m=2, t0=0.001,tend=0.01, unknown_storage=:dense,solver=nothing)
    sys,X=create_porous_medium_problem(n,m,unknown_storage)
    inival=unknowns(sys)
    inival[1,:].=map(x->barenblatt(x,t0,m),X)
    problem = ODEProblem(sys,inival,(t0,tend))
    odesol = DifferentialEquations.solve(problem,Rodas5(linsolve=UMFPACKFactorization()))
    sol=reshape(odesol,sys)
    err=norm(sol[1,:,end]-map(x->barenblatt(x,tend,m),X))
    sol, sys,err
end


function main(;m=2,n=20, solver=nothing, unknown_storage=:dense, Plotter=nothing)

    vis=GridVisualizer(Plotter=Plotter,layout=(1,2),resolution=(800,400))
    t=@elapsed begin
        sol1,sys,err1=run_vfvm(m=m,n=n, unknown_storage=unknown_storage)
    end
    println(history_summary(sys))
    title=@sprintf("VoronoiFVM: %.0f ms e=%.2e",t*1000,err1)
    println(title)
    scalarplot!(vis[1,1],sys,sol1,title=title,aspect=400)

    t=@elapsed begin
        sol2,sys,err2=run_diffeq(m=m,n=n,solver=solver, unknown_storage=unknown_storage)
    end
    println(history_summary(sys))
    title=@sprintf("    DiffEq: %.0f ms, e=%.2e",t*1000,err2)
    println(title)
    scalarplot!(vis[1,2],sys,sol2,title=title,aspect=400)
    reveal(vis)
    norm(sol2[end]-sol1[end],Inf)<0.01

end

#test()=main()

test()=true

end

