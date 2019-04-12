module ThreeRegions1D

using Printf
using VoronoiFVM

if isinteractive()
    using PyPlot
end




function main(;n=30,pyplot=false,verbose=false,dense=false)
    h=3.0/(n-1)
    X=collect(0:h:3.0)
    grid=VoronoiFVM.Grid(X)
    cellmask!(grid,[0.0],[1.0],1)
    cellmask!(grid,[1.0],[2.1],2)
    cellmask!(grid,[1.9],[3.0],3)

    subgrid1=subgrid(grid,[1])   
    subgrid2=subgrid(grid,[1,2,3])
    subgrid3=subgrid(grid,[3])
    
    if pyplot
        clf()
        fvmplot(grid)
        show()
        waitforbuttonpress()

    end
    
    eps=[1,1,1]
    k=[1,1,1]

    physics=VoronoiFVM.Physics(
    num_species=3,
    reaction=function(f,u,node,data)
        if node.region==1
            f[1]=k[1]*u[1]
            f[2]=-k[1]*u[1]
        elseif node.region==3
            f[2]=k[3]*u[2]
            f[3]=-k[3]*u[2]
        else
            f[1]=0
        end
    end,
    
    flux=function(f,u,edge,data)   
        uk=viewK(edge,u)
        ul=viewL(edge,u)
        if edge.region==1
            f[1]=eps[1]*(uk[1]-ul[1])
            f[2]=eps[2]*(uk[2]-ul[2])
        elseif edge.region==2
            f[2]=eps[2]*(uk[2]-ul[2])
        elseif edge.region==3
            f[2]=eps[2]*(uk[2]-ul[2])
            f[3]=eps[3]*(uk[3]-ul[3])
        end
    end,
    
    source=function(f,node,data)
        if node.region==1
            f[1]=1.0e-4*(3.0-node.coord[1])
        end
    end,
    
    storage=function(f,u,node,data)
        f.=u
    end
    )

    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics)
    end

    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1,2,3])
    enable_species!(sys,3,[3])

    sys.boundary_factors[3,2]=VoronoiFVM.Dirichlet
    sys.boundary_values[3,2]=0
    
    inival=unknowns(sys)
    U=unknowns(sys)
    inival.=0
    
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    tstep=0.01
    time=0.0
    istep=0
    tend=10
    testval=0
    while time<tend
        time=time+tstep
        solve!(U,inival,sys,control=control,tstep=tstep)
        inival.=U
        if verbose
            @printf("time=%g\n",time)
        end
        tstep*=1.0
        istep=istep+1
        testval=U[2,15]
        if pyplot && istep%10 == 0
            PyPlot.clf()
            fvmplot(subgrid1, U[1,:],label="spec1", color=(0.5,0,0))
            fvmplot(subgrid2, U[2,:],label="spec2", color=(0.0,0.5,0))
            fvmplot(subgrid3, U[3,:],label="spec3", color=(0.0,0.0,0.5))
            PyPlot.legend(loc="best")
            PyPlot.grid()
            pause(1.0e-10)
        end
    end
    return testval
end
end

