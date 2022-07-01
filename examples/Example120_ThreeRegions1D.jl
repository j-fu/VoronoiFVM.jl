# # 120: Differing species sets in regions, 1D
# ([source code](SOURCE_URL))

module Example120_ThreeRegions1D

using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize



function main(;n=30,Plotter=nothing,plot_grid=false, verbose=false,unknown_storage=:sparse,tend=10)
    h=3.0/(n-1)
    X=collect(0:h:3.0)
    grid=VoronoiFVM.Grid(X)
    cellmask!(grid,[0.0],[1.0],1)
    cellmask!(grid,[1.0],[2.1],2)
    cellmask!(grid,[1.9],[3.0],3)

    subgrid1=subgrid(grid,[1])   
    subgrid2=subgrid(grid,[1,2,3])
    subgrid3=subgrid(grid,[3])
    
    if plot_grid
        plotgrid(grid,Plotter=Plotter)
        return
    end
    
    eps=[1,1,1]
    k=[1,1,1]

    physics=VoronoiFVM.Physics(
    reaction=function(f,u,node)
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
    
        flux=function(f,u,edge)
        
        if edge.region==1
            f[1]=eps[1]*(u[1,1]-u[1,2])
            f[2]=eps[2]*(u[2,1]-u[2,2])
        elseif edge.region==2
            f[2]=eps[2]*(u[2,1]-u[2,2])
        elseif edge.region==3
            f[2]=eps[2]*(u[2,1]-u[2,2])
            f[3]=eps[3]*(u[3,1]-u[3,2])
        end
    end,
    
    source=function(f,node)
        if node.region==1
           f[1]=1.0e-4*(3.0-node[1])
        end
    end,
    
    storage=function(f,u,node)
       if node.region==1
            f[1]=u[1]
            f[2]=u[2]
        elseif node.region==2
            f[2]=u[2]
        elseif node.region==3
            f[2]=u[2]
            f[3]=u[3]
        end
    end
    )

    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)

    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1,2,3])
    enable_species!(sys,3,[3])

    boundary_dirichlet!(sys,3,2,0.0)
    
    inival=unknowns(sys)
    U=unknowns(sys)
    inival.=0
    
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    tstep=0.01
    time=0.0
    istep=0
    testval=0
    p=GridVisualizer(Plotter=Plotter,layout=(1,1))
    while time<tend
        time=time+tstep
        solve!(U,inival,sys,control=control,tstep=tstep)
        inival.=U
        if verbose
            @printf("time=%g\n",time)
        end
        tstep*=1.1
        istep=istep+1
        testval=U[2,5]
        
        U1=view(U[1,:],subgrid1)
        U2=view(U[2,:],subgrid2)
        U3=view(U[3,:],subgrid3)
        
        scalarplot!(p[1,1],subgrid1, U1,label="spec1", color=(0.5,0,0),xlimits=(0,3),flimits=(0,1e-3),
                  title=@sprintf("three regions t=%.3g",time))
        scalarplot!(p[1,1],subgrid2, U2,label="spec2", color=(0.0,0.5,0),clear=false)
        scalarplot!(p[1,1],subgrid3, U3,label="spec3", color=(0.0,0.0,0.5),clear=false,show=true)
        sleep(1.0e-2)
    end
    return testval
end

function test()
    testval=0.0005967243505359461
    main(unknown_storage=:sparse) ≈ testval &&
        main(unknown_storage=:dense) ≈ testval
end

end
