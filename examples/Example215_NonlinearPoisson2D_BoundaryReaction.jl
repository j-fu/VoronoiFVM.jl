# # 215: 2D Nonlinear Poisson with boundary reaction

module Example215_NonlinearPoisson2D_BoundaryReaction

using Printf
using VoronoiFVM
const Node=VoronoiFVM.Node
const Edge=VoronoiFVM.Edge

if installed("Plots")
    using Plots
end

function main(;n=10,doplot=false,verbose=false, dense=false)
    if !installed("Plots")
        doplot=false
    end
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)
    
    grid=VoronoiFVM.Grid(X,Y)

    
    eps=1.0e-2
    k=1.0
    physics=VoronoiFVM.Physics(
        num_species=2,
        breaction=function(f,u,node,data)
        if  node.region==2
            f[1]=k*(u[1]-u[2])
            f[2]=k*(u[2]-u[1])
        else
            f[1]=0        
            f[2]=0
        end
        end,
    
    flux=function(f,u,edge,data)
        uk=viewK(2,u)
        ul=viewL(2,u)
        f[1]=eps*(uk[1]-ul[1])
        f[2]=eps*(uk[2]-ul[2])
    end,
    
    source=function(f,node,data)
        x1=node.coord[1]-0.5
        x2=node.coord[2]-0.5
        f[1]=exp(-20.0*(x1^2+x2^2))
    end,
    
    storage=function(f,u,node,data)
        f[1]=u[1]
        f[2]=u[2]
    end)

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
        @views if doplot
            p1=contourf(X,Y,reshape(U[1,:],length(X),length(Y)),levels=collect(0:0.125:0.75),clim=(0,0.75),colorbar=:right,color=:viridis,title=@sprintf("max1=%g max2=%g\n",maximum(U[1,:]),maximum(U[2,:])))
            p2=contourf(X,Y,reshape(U[2,:],length(X),length(Y)),levels=collect(0:0.0025:0.02),clim=(0,0.02), colorbar=:right,color=:viridis)
            p=Plots.plot(p1,p2,layout=(2,1) )
            gui(p)
        end
    end
    return u25
end

function test()
    main() ≈ 0.008761335823958986 &&
    main(dense=true) ≈ 0.008761335823958986
end
end
