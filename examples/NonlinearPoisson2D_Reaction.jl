module NonlinearPoisson2D_Reaction

using Printf
using TwoPointFluxFVM

if isinteractive()
    using PyPlot
end

mutable struct Physics <: TwoPointFluxFVM.Physics
    reaction::Function
    flux::Function
    source::Function
    storage::Function
    eps::Float64
    k::Float64
    Physics()=new()
end


function main(;n=10,pyplot=false,verbose=false, dense=false)
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)

    grid=TwoPointFluxFVM.Grid(X,Y)
    
    
    physics=Physics()
    physics.eps=1.0e-2
    physics.k=1.0
    
    physics.reaction=function(physics,node,f,u)
        f[1]=physics.k*(u[1]-u[2])
        f[2]=physics.k*(u[2]-u[1])
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
        sys=TwoPointFluxFVM.DenseSystem(grid,physics,2)
    else
        sys=TwoPointFluxFVM.SparseSystem(grid,physics,2)
    end

    add_species(sys,1,[1])
    add_species(sys,2,[1])
    
    inival=unknowns(sys)
    inival.=0.0
    
    
    control=TwoPointFluxFVM.NewtonControl()
    control.verbose=verbose
    control.tol_linear=1.0e-5
    control.max_lureuse=0
    tstep=0.01
    time=0.0
    istep=0
    u15=0
    while time<1
        time=time+tstep
        U=solve(sys,inival,control=control,tstep=tstep)
        inival.=U
        if verbose
            @printf("time=%g\n",time)
        end
        u15=U[15]
        tstep*=1.0
        istep=istep+1
        
        @views if pyplot && istep%10==1
            if verbose
                @printf("max1=%g max2=%g\n",maximum(U[1,:]),maximum(U[2,:]))
            end
            levels1=collect(0:0.01:1)
            PyPlot.clf()
            subplot(211)
            contourf(X,Y,reshape(U[1,:],length(X),length(Y)), cmap=ColorMap("hot"),vmin=0, vmax=0.6)
            colorbar()
            levels2=collect(0:0.001:0.1)
            subplot(212)
            contourf(X,Y,reshape(U[2,:],length(X),length(Y)), cmap=ColorMap("hot"),vmin=0, vmax=0.6)
            colorbar()
            
            pause(1.0e-10)
        end
    end
    return u15
end

end
