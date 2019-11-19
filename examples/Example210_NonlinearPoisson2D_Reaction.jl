# # 210: 2D Nonlinear Poisson with reaction

module Example210_NonlinearPoisson2D_Reaction

using Printf
using VoronoiFVM

if installed("Plots")
    using Plots
end




mutable struct MyData <: VoronoiFVM.AbstractData
    eps::Float64
    k::Float64
    MyData()=new() 
end



function main(;n=10,doplot=false,verbose=false, dense=false)
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)

    grid=VoronoiFVM.Grid(X,Y)
    data=MyData()
    
    function reaction!(f,u,node,data)
        f[1]=data.k*(u[1]-u[2])
        f[2]=data.k*(u[2]-u[1])
    end
    
    function flux!(f,u,edge,data)
        uk=viewK(2,u)
        ul=viewL(2,u)
        f[1]=data.eps*(uk[1]-ul[1])
        f[2]=data.eps*(uk[2]-ul[2])
    end
    
    function source!(f,node,data)
        x1=node.coord[1]-0.5
        x2=node.coord[2]-0.5
        f[1]=exp(-20*(x1^2+x2^2))
    end
    
    function storage!(f,u,node,data)
        f[1]=u[1]
        f[2]=u[2]
    end
    
    
    physics=VoronoiFVM.Physics(num_species=2,
                               data=data,
                               flux=flux!,
                               storage=storage!,
                               reaction=reaction!,
                               source=source!)
    
    data.eps=1.0e-2
    data.k=1.0

    
    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics)
    end

    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])
    
    inival=unknowns(sys)
    U=unknowns(sys)
    inival.=0.0
    
    
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    control.tol_linear=1.0e-5
    control.max_lureuse=0
    tstep=0.01
    time=0.0
    istep=0
    u15=0
    while time<1
        time=time+tstep
        solve!(U,inival,sys,control=control,tstep=tstep)
        inival.=U
        if verbose
            @printf("time=%g\n",time)
        end
        u15=U[15]
        tstep*=1.0
        istep=istep+1
        
        @views if doplot
            p1=contourf(X,Y,reshape(U[1,:],length(X),length(Y)),levels=collect(0:0.1:0.6),clim=(0,0.6),colorbar=:right,color=:viridis,title=@sprintf("max1=%g max2=%g\n",maximum(U[1,:]),maximum(U[2,:])))
            p2=contourf(X,Y,reshape(U[2,:],length(X),length(Y)),levels=collect(0:0.1:0.6),clim=(0,0.6), colorbar=:right,color=:viridis)
            p=Plots.plot(p1,p2,layout=(2,1) )
            gui(p)
        end
    end
    return u15
end

function test()
    main(dense=false) ≈ 0.014566189535134827 &&
        main(dense=true) ≈ 0.014566189535134827
end
end
