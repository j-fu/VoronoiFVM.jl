# 115: 1D two species system with boundary reaction and boundary species


module Example115_NonlinearPoisson1D_BoundarySpecies

using Printf
using VoronoiFVM
const Node=VoronoiFVM.Node
const Edge=VoronoiFVM.Edge

if isinteractive()
    using Plots
end



function main(;n=10,doplot=false,verbose=false,tend=1, dense=false)
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    N=length(X)
    
    grid=VoronoiFVM.Grid(X)

    
    eps=21
    k=1
    
    physics=VoronoiFVM.Physics(
        num_species=3,
        breaction=function(f,u,node,data)
        if  node.region==2
        f[1]=k*(u[1]-u[3])
        f[2]=k*(u[2]-u[3])
        f[3]=k*(u[3]-u[1])+ k*(u[3]-u[2])
        end
        end,
        
        bstorage=function(f,u,node,data)
        if  node.region==2
        f[3]=u[3]
        end
        end,
        
        flux=function(f,u,edge,data)
        uk=viewK(edge,u)
        ul=viewL(edge,u)
        f[1]=eps*(uk[1]-ul[1])
        f[2]=eps*(uk[2]-ul[2])
        end ,
        
        source=function(f,node,data)
        x1=node.coord[1]-0.5
        f[1]=exp(-20*x1^2)
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
    enable_boundary_species!(sys,3,[2])

    inival=unknowns(sys)
    inival.=0.0
    U=unknowns(sys)
    
    eps=1.0e-2
    
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    control.tol_linear=1.0e-5
    control.tol_relative=1.0e-5
    control.max_lureuse=0
    tstep=0.01
    time=0.0
    istep=0
    T=zeros(0)
    Ub=zeros(0)
    u5=0
    while time<tend
        time=time+tstep
        solve!(U,inival,sys, control=control,tstep=tstep)
        inival.=U
        if verbose
            @printf("time=%g\n",time)
        end
        tstep*=1.0
        istep=istep+1
        u5=U[1,3]
        
        append!(T,time)
        append!(Ub,U[3,N])
        
        if doplot
            p1=Plots.plot(grid.coord[1,:],U[1,:], grid=true, label="spec1")
            Plots.plot!(p1,grid.coord[1,:],U[2,:], label="spec2",title=@sprintf("max1=%.5f max2=%.5f maxb=%.5f\n",maximum(U[1,:]),maximum(U[2,:]),U[3,N]))
            p2=Plots.plot(T,Ub,ylabel="U_b",xlabel="t")
            p=Plots.plot(p1,p2,layout=(2,1),legend=false)
            gui(p)
        end
    end
    return u5
end

function test()
    main(dense=false) ≈ 0.22631106953924143 &&
        main(dense=true) ≈ 0.22631106953924143
end
end
