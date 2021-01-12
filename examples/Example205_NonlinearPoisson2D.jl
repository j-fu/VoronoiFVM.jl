# # 205: 2D Nonlinear Poisson equation
# ([source code](SOURCE_URL))

module Example205_NonlinearPoisson2D

using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize



function main(;n=10,Plotter=nothing,verbose=false, unknown_storage=:sparse)
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)


    grid=VoronoiFVM.Grid(X,Y)
    
    eps=1.0e-2
    
    physics=VoronoiFVM.Physics(
        num_species=1,
        reaction=function(f,u,node)
        f[1]=u[1]^2
        end,
        
        flux=function(f,u0,edge)
        u=unknowns(edge,u0)
        f[1]=eps*(u[1,1]^2-u[1,2]^2)
        end,
        
        source=function(f,node)
        x1=node[1]-0.5
        x2=node[2]-0.5
        f[1]=exp(-20.0*(x1^2+x2^2))
        end,
        
        storage=function(f,u,node)
        f[1]=u[1]
        end)
    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)
    enable_species!(sys,1,[1])

    boundary_dirichlet!(sys,1,2,0.1)
    boundary_dirichlet!(sys,1,4,0.1)
    
    inival=unknowns(sys)
    U=unknowns(sys)
    inival.=0.5


    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    control.tol_linear=1.0e-5
    control.max_lureuse=10
    tstep=0.01
    time=0.0
    u15=0
    p=GridVisualizer(Plotter=Plotter)
    while time<1.0
        time=time+tstep
        solve!(U,inival,sys,control=control,tstep=tstep)
        u15=U[15]
        inival.=U

        if verbose
            @printf("time=%g\n",time)
        end

        scalarplot!(p[1,1],grid,U[1,:],Plotter=Plotter,clear=true,show=true)
        tstep*=1.0
    end
    return u15
end

function test()
    main(unknown_storage=:sparse) ≈ 0.3554284760906605 &&
        main(unknown_storage=:dense) ≈ 0.3554284760906605
end
end
