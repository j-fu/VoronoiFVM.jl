module NonlinearPoisson2D

using Printf
using VoronoiFVM

if isinteractive()
    using PyPlot
end



function main(;n=10,pyplot=false,verbose=false, dense=false)
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)


    grid=VoronoiFVM.Grid(X,Y)
    
    eps=1.0e-2
    
    physics=VoronoiFVM.Physics(
        num_species=1,
        reaction=function(f,u,node,data)
        f[1]=u[1]^2
        end,
        
        flux=function(f,u,edge,data)
        f[1]=eps*(u[1]^2-u[2]^2)
        end,
        
        source=function(f,node,data)
        x1=node.coord[1]-0.5
        x2=node.coord[2]-0.5
        f[1]=exp(-20.0*(x1^2+x2^2))
        end,
        
        storage=function(f,u,node,data)
        f[1]=u[1]
        end)
    
    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics)
    end        
    enable_species!(sys,1,[1])

    sys.boundary_values[1,2]=0.1
    sys.boundary_values[1,4]=0.1
    
    sys.boundary_factors[1,2]=VoronoiFVM.Dirichlet
    sys.boundary_factors[1,4]=VoronoiFVM.Dirichlet
    
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
    while time<1.0
        time=time+tstep
        solve!(U,inival,sys,control=control,tstep=tstep)
        u15=U[15]
        inival.=U

        if verbose
            @printf("time=%g\n",time)
        end

        tstep*=1.0
        if pyplot
            levels=collect(0:0.01:1)
            PyPlot.clf()
            contourf(X,Y,reshape(values(U),length(X),length(Y)), cmap=ColorMap("hot"),levels=levels)
            colorbar()
            pause(1.0e-10)
        end
    end
    return u15
end

end
