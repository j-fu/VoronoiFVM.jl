module NonlinearPoisson2D_BoundarySpecies

using Printf
using TwoPointFluxFVM


if isinteractive()
    using PyPlot
end


mutable struct Physics <: TwoPointFluxFVM.Physics
    flux::Function
    source::Function
    storage::Function
    bstorage::Function
    breaction::Function
    k::Float64
    eps::Float64 
    Physics()=new()
end



function main(;n=10,pyplot=false,verbose=false,dense=false)
    
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)
    


    grid=TwoPointFluxFVM.Grid(X,Y)
    
    
    physics=Physics()
    physics.k=1
    physics.eps=1
    
    physics.breaction=function(physics,node,f,u)
        if  node.region==2
            f[1]=physics.k*(u[1]-u[3])
            f[3]=physics.k*(u[3]-u[1])+ physics.k*(u[3]-u[2])
            f[2]=physics.k*(u[2]-u[3])
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
        x2=node.coord[2]-0.5
        f[1]=exp(-20*(x1^2+x2^2))
    end 
    
    physics.storage=function(physics,node, f,u)
        f[1]=u[1]
        f[2]=u[2]
    end
    
    if dense
        sys=TwoPointFluxFVM.SparseSystem(grid,physics,3)
    else
        sys=TwoPointFluxFVM.SparseSystem(grid,physics,3)
    end
    add_species(sys,1,[1])
    add_species(sys,2,[1])
    add_boundary_species(sys,3,[2])

    
    function tran32!(a,b)
        a[1]=b[2]
    end
    
    bgrid2=subgrid(grid,[2],boundary=true,transform=tran32!)
   
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
    u5=0
    while time<1
        time=time+tstep
        U=solve(sys,inival,control=control,tstep=tstep)
        inival.=U
        if verbose
            @printf("time=%g\n",time)
        end
        tstep*=1.0
        istep=istep+1
        U_bound=view(U,bgrid2)
        u5=U_bound[3,5]
        if pyplot && istep%10 == 0
            @printf("max1=%g max2=%g maxb=%g\n",maximum(U[1,:]),maximum(U[2,:]),maximum(U_bound[3,:]))
            PyPlot.clf()

            subplot(311)
            levels1=collect(0:0.01:1)
            contourf(X,Y,reshape(U[1,:],length(X),length(Y)), cmap=ColorMap("hot"))
            colorbar()

            subplot(312)
            levels2=collect(0:0.001:0.1)
            contourf(X,Y,reshape(U[2,:],length(X),length(Y)), cmap=ColorMap("hot"))
            colorbar()

            subplot(313)
            levels2=collect(0:0.001:0.1)
            fvmplot(bgrid2,U[3,:])
            pause(1.0e-10)
        end
    end
    return u5
end
end
