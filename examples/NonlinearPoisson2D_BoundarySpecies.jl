module NonlinearPoisson2D_BoundarySpecies

using Printf
using TwoPointFluxFVM
if isinteractive()
    using PyPlot
end


mutable struct Physics <:TwoPointFluxFVM.Physics
    TwoPointFluxFVM.@AddPhysicsBaseClassFields
    k::Float64
    eps::Float64 
    Physics()=Physics(new())
end


function breaction!(this::Physics,f,bf,u,bu)
    if  this.bregion==2
        f[1]=this.k*(u[1]-bu[1])
        bf[1]=this.k*(bu[1]-u[1])+ this.k*(bu[1]-u[2])
        f[2]=this.k*(u[2]-bu[1])
    else
        f[1]=0        
        f[2]=0
    end
end

function bstorage!(this::Physics,bf,bu)
    if  this.bregion==2
        bf[1]=bu[1]
    else
        bf[1]=0        
    end
end


function flux!(this::Physics,f,uk,ul)
    f[1]=this.eps*(uk[1]-ul[1])
    f[2]=this.eps*(uk[2]-ul[2])
end 

function source!(this::Physics,f,x)
    x1=x[1]-0.5
    x2=x[2]-0.5
    f[1]=exp(-20*(x1^2+x2^2))
end 



function Physics(this)
    TwoPointFluxFVM.PhysicsBase(this,2)
    this.num_bspecies=[ 0, 1, 0, 0]
    this.eps=1
    this.k=1.0
    this.breaction=breaction!
    this.bstorage=bstorage!
    this.flux=flux!
    this.source=source!
    return this
end


function main(;n=10,pyplot=false,verbose=false)
    
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)
    


    geom=TwoPointFluxFVM.Graph(X,Y)

    
    physics=Physics()
    
    
    sys=TwoPointFluxFVM.System(geom,physics)
    
    inival=unknowns(sys)
    inival.=0.0
    
    physics.eps=1.0e-2
    
    control=TwoPointFluxFVM.NewtonControl()
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
        U=solve(sys,inival,control=control,tstep=tstep)
        inival.=U
        # for i in eachindex(U)
        #     inival[i]=U[i]
        # end
        if verbose
            @printf("time=%g\n",time)
        end
        tstep*=1.0
        istep=istep+1
        U_bulk=bulk_unknowns(sys,U)
        U_bound=boundary_unknowns(sys,U,2)
        u5=U_bound[5]
        if pyplot && istep%10 == 0
            @printf("max1=%g max2=%g maxb=%g\n",maximum(U_bulk[1,:]),maximum(U_bulk[2,:]),maximum(U_bound))
            PyPlot.clf()

            subplot(311)
            levels1=collect(0:0.01:1)
            contourf(X,Y,reshape(U_bulk[1,:],length(X),length(Y)), cmap=ColorMap("hot"))
            colorbar()

            subplot(312)
            levels2=collect(0:0.001:0.1)
            contourf(X,Y,reshape(U_bulk[2,:],length(X),length(Y)), cmap=ColorMap("hot"))
            colorbar()

            subplot(313)
            levels2=collect(0:0.001:0.1)
            plot(X,U_bound[1,:])
            pause(1.0e-10)
        end
    end
    return u5
end
end
