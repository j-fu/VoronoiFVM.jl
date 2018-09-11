using Printf
using TwoPointFluxFVM

if !isinteractive()
    using PyPlot
end


mutable struct ILiqParameters <:FVMParameters
    number_of_species::Int64
    eps::Float64 
    a::Float64
    function ILiqParameters()
        new(2,1.0e-4,0)
    end
end

const iphi=1
const ic=2
const beps=1.0e-4

function run_iliq(;n=100,pyplot=false)

    h=1.0/convert(Float64,n)
    geom=FVMGraph(collect(0:h:1))
    
    parameters=ILiqParameters()
    
    function flux!(this::ILiqParameters,f,uk,ul)
        f[iphi]=this.eps*(uk[iphi]-ul[iphi])
        muk=-log(1-uk[ic])
        mul=-log(1-ul[ic])
        bp,bm=fbernoulli_pm(2*(uk[iphi]-ul[iphi])+(muk-mul))
        f[ic]=bm*uk[ic]-bp*ul[ic]
    end 


    function classflux!(this::ILiqParameters,f,uk,ul)
        f[iphi]=this.eps*(uk[iphi]-ul[iphi])
        arg=uk[iphi]-ul[iphi]
        bp,bm=fbernoulli_pm(uk[iphi]-ul[iphi])
        f[ic]=bm*uk[ic]-bp*ul[ic]
    end 

    function storage!(this::FVMParameters, f,u)
        f[iphi]=0
        f[ic]=u[ic]
    end
    
    function reaction!(this::FVMParameters, f,u)
        f[iphi]=1-2*u[ic]
        f[ic]=0
    end
    
    
    sys=TwoPointFluxFVMSystem(geom,parameters=parameters, 
                              storage=storage!, 
                              flux=flux!, 
                              reaction=reaction!
                              )
    sys.boundary_values[iphi,1]=1
    sys.boundary_values[iphi,2]=0.0
    
    sys.boundary_factors[iphi,1]=Dirichlet
    sys.boundary_factors[iphi,2]=Dirichlet
    
    inival=unknowns(sys)
    for inode=1:size(inival,2)
        inival[iphi,inode]=0
        inival[ic,inode]=0.6
    end
    parameters.eps=1.0e-2
    parameters.a=5
    control=FVMNewtonControl()
    control.verbose=true
    t=0.0
    tend=1.0
    tstep=1.0e-10
    while t<tend
        t=t+tstep
        U=solve(sys,inival,control=control,tstep=tstep)
        for i=1:size(inival,2)
            inival[iphi,i]=U[iphi,i]
            inival[ic,i]=U[ic,i]
        end
        @printf("time=%g\n",t)
        if pyplot
            PyPlot.clf()
            plot(geom.Nodes[1,:],U[iphi,:])
            plot(geom.Nodes[1,:],U[ic,:])
            pause(1.0e-10)
        end
        tstep*=1.2
    end
end



if !isinteractive()
    @time run_iliq(n=100,pyplot=true)
    waitforbuttonpress()
end
