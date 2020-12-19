# # 215: 2D Nonlinear Poisson with boundary reaction
# ([source code](SOURCE_URL))

module Example215_NonlinearPoisson2D_BoundaryReaction

using Printf
using VoronoiFVM

function main(;n=10,Plotter=nothing,verbose=false, unknown_storage=:sparse)
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)
    
    grid=VoronoiFVM.Grid(X,Y)

    
    eps=1.0e-2
    k=1.0
    physics=VoronoiFVM.Physics(
        num_species=2,
        breaction=function(f,u,node)
        if  node.region==2
            f[1]=k*(u[1]-u[2])
            f[2]=k*(u[2]-u[1])
        else
            f[1]=0        
            f[2]=0
        end
        end,
    
    flux=function(f,u0,edge)
        u=unknowns(edge,u0)
        f[1]=eps*(u[1,1]-u[1,2])
        f[2]=eps*(u[2,1]-u[2,2])
    end,
    
    source=function(f,node)
        x1=node[1]-0.5
        x2=node[2]-0.5
        f[1]=exp(-20.0*(x1^2+x2^2))
    end,
    
    storage=function(f,u,node)
        f[1]=u[1]
        f[2]=u[2]
    end)

    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)
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
    u25=0
    while time<1
        time=time+tstep
        solve!(U,inival,sys,control=control,tstep=tstep)
        inival.=U
        # for i in eachindex(U)
        #     inival[i]=U[i]
        # end
        if verbose
            @printf("time=%g\n",time)
        end

        tstep*=1.0
        istep=istep+1
        u25=U[25]
        if isplots(Plotter)
            p1=Plotter.contourf(X,Y,reshape(U[1,:],length(X),length(Y)),levels=collect(0:0.125:0.75),clim=(0,0.75),colorbar=:right,color=:viridis,title=@sprintf("max1=%g max2=%g\n",maximum(U[1,:]),maximum(U[2,:])))
            p2=Plotter.contourf(X,Y,reshape(U[2,:],length(X),length(Y)),levels=collect(0:0.0025:0.02),clim=(0,0.02), colorbar=:right,color=:viridis)
            p=Plotter.plot(p1,p2,layout=(2,1) )
            Plotter.gui(p)
        end
    end
    return u25
end

function test()
    main(unknown_storage=:sparse) ≈ 0.008761335823958986 &&
    main(unknown_storage=:dense) ≈ 0.008761335823958986
end
end
