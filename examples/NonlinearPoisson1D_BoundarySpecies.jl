module NonlinearPoisson1D_BoundarySpecies

using Printf
using TwoPointFluxFVM
const Node=TwoPointFluxFVM.Node
const Edge=TwoPointFluxFVM.Edge

if isinteractive()
    using PyPlot
end


mutable struct Physics
    flux::Function
    source::Function
    storage::Function
    bstorage::Function
    breaction::Function
    
    k::Float64
    eps::Float64 
    Physics()=new()
end



function main(;n=10,pyplot=false,verbose=false,tend=1)
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    N=length(X)
    
    grid=TwoPointFluxFVM.Grid(X)

    
    physics=Physics()
    physics.eps=21
    physics.k=1
    

    physics.breaction=function(physics,node,f,u)
        if  node.region==2
            f[1]=physics.k*(u[1]-u[3])
            f[2]=physics.k*(u[2]-u[3])
            f[3]=physics.k*(u[3]-u[1])+ physics.k*(u[3]-u[2])
        end
    end
    
    physics.bstorage=function(physics,node,f,u)
        if  node.region==2
            f[3]=u[3]
        end
    end
    
    physics.flux=function(physics,edge,f,uk,ul)
        f[1]=physics.eps*(uk[1]-ul[1])
        f[2]=physics.eps*(uk[2]-ul[2])
    end 
    
    physics.source=function(physics,node,f)
        x1=node.coord[1]-0.5
        f[1]=exp(-20*x1^2)
    end 
    physics.storage=function(physics,node, f,u)
        f[1]=u[1]
        f[2]=u[2]
    end
    

    sys=TwoPointFluxFVM.System(grid,physics,3)
    add_species(sys,1,[1])
    add_species(sys,2,[1])
    add_boundary_species(sys,3,[2])

    inival=unknowns(sys)
    inival.=0.0
    
    physics.eps=1.0e-2
    
    control=TwoPointFluxFVM.NewtonControl()
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
        U=solve(sys,inival,control=control,tstep=tstep)
        inival.=U
        if verbose
            @printf("time=%g\n",time)
        end
        tstep*=1.0
        istep=istep+1
        u5=getdof(U,5)
        
        append!(T,time)
        append!(Ub,U[3,N])
        
        if pyplot && istep%10 == 0
            @printf("max1=%g max2=%g maxb=%g\n",maximum(U[1,:]),maximum(U[2,:]),U[3,N])
            PyPlot.clf()
            subplot(211)
            plot(X,U[1,:],label="spec1")
            plot(X,U[2,:],label="spec2")
            PyPlot.legend(loc="best")
            PyPlot.grid()
            subplot(212)
            plot(T,Ub,label="U_b")
            PyPlot.legend(loc="best")
            PyPlot.grid()
            pause(1.0e-10)
        end
    end
    return u5
end
end
