# # 220: 2D Nonlinear Poisson with boundary reaction and boundary species
# ([source code](SOURCE_URL))

module Example220_NonlinearPoisson2D_BoundarySpecies

using Printf
using VoronoiFVM


function main(;n=10,Plotter=nothing,verbose=false,dense=false)
    
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)
    


    grid=VoronoiFVM.Grid(X,Y)
    
    
    k=1.0
    eps=1.0
    physics=VoronoiFVM.Physics(
    num_species=3,
    breaction=function(f,u,node)
        if  node.region==2
            f[1]=k*(u[1]-u[3])
            f[3]=k*(u[3]-u[1])+ k*(u[3]-u[2])
            f[2]=k*(u[2]-u[3])
        end
    end,
    
    bstorage=function(f,u,node)
        if  node.region==2
            f[3]=u[3]
        end
    end,
    
    
    flux=function(f,u,edge)
        uk=viewK(edge,u)
        ul=viewL(edge,u)
        f[1]=eps*(uk[1]-ul[1])
        f[2]=eps*(uk[2]-ul[2])
    end,
    
    source=function(f,node)
        x1=node[1]-0.5
        x2=node[2]-0.5
        f[1]=exp(-20.0*(x1^2+x2^2))
    end,
    
    storage=function(f,u,node)
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

    
    function tran32!(a,b)
        a[1]=b[2]
    end
    
    bgrid2=subgrid(grid,[2],boundary=true,transform=tran32!)
   
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
    u5=0
    while time<1
        time=time+tstep
        solve!(U,inival,sys,control=control,tstep=tstep)
        inival.=U
        if verbose
            @printf("time=%g\n",time)
        end
        tstep*=1.0
        istep=istep+1
        U_bound=view(U,bgrid2)
        u5=U_bound[3,5]
        if isplots(Plotter)
            p1=Plotter.contourf(X,Y,reshape(U[1,:],length(X),length(Y)),levels=collect(0:0.1:0.6),clim=(0,0.6),colorbar=:right,color=:viridis,title=@sprintf("max1=%g max2=%g maxb=%g\n",maximum(U[1,:]),maximum(U[2,:]),maximum(U_bound[3,:])))
            p2=Plotter.contourf(X,Y,reshape(U[2,:],length(X),length(Y)),levels=collect(0:0.0001:0.002),clim=(0,0.002), colorbar=:right,color=:viridis)
            p3=Plotter.plot(grid=true,ylims=(0,0.0025))
            VoronoiFVM.plot(Plotter,bgrid2,U[3,:],p=p3,show=false)
            p=Plotter.plot(p1,p2,p3,layout=(3,1) )
            Plotter.gui(p)
        end
    end
    return u5
end

function test()
    main() ≈ 0.0020781361856598
    main(dense=true) ≈ 0.0020781361856598
end
end
