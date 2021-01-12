# # 210: 2D Nonlinear Poisson with reaction
# ([source code](SOURCE_URL))

module Example210_NonlinearPoisson2D_Reaction

using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize

function main(;n=10,Plotter=nothing,verbose=false, unknown_storage=:sparse)

    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)

    grid=VoronoiFVM.Grid(X,Y)
    data = (eps=1.0e-2, k=1.0)

    function reaction!(f,u,node,data)
        f[1]=data.k*(u[1]-u[2])
        f[2]=data.k*(u[2]-u[1])
    end

    function flux!(f,u0,edge,data)
        u=unknowns(edge,u0)
        f[1]=data.eps*(u[1,1]-u[1,2])
        f[2]=data.eps*(u[2,1]-u[2,2])
    end

    function source!(f,node,data)
        x1=node[1]-0.5
        x2=node[2]-0.5
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
    u15=0
    p=GridVisualizer(Plotter=Plotter,layout=(2,1))
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
        scalarplot!(p[1,1],grid,U[1,:],clear=true)
        scalarplot!(p[2,1],grid,U[2,:],show=true)
    end
    return u15
end

function test()
    main(unknown_storage=:sparse) ≈ 0.014566189535134827 &&
        main(unknown_storage=:dense) ≈ 0.014566189535134827
end
end
