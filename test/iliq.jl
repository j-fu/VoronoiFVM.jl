using Printf
if !isinteractive()
    # Command line argument parser, get it via Pkg.add
    using TwoPointFluxFVM
    using ArgParse
    argdef = ArgParseSettings()
    add_arg_table(argdef,
                  "--pyplot", Dict(:help => "call python visualization",
                                   #                    :default => false,
                                   :action => :store_true),
                  "--n", Dict(:help => "problem size",
                              :arg_type => Int,
                              :default => 11)
                  
                  )
    
    args = parse_args(argdef)
    
    
    pyplot=false
    if args["pyplot"]
        using PyPlot
        pyplot=true
    end
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

function run_iliq(n;pyplot=false)
    
    
    
    function bernoulli(x)
        if x<-beps
            return x/(exp(x)-1)
        elseif x <beps
            x2  = x*x;
            x4  = x2*x2;
            x6  = x4*x2;
            x8  = x6*x2;
            x10 = x8*x2;
            return 1.0 - 0.5*x +1.0/12.0 * x2 
            - 1.0/720.0 * x4 
            + 1.0/30240.0 * x6 
            - 1.0/1209600.0 * x8 
            + 1.0/47900160.0 * x10; 
        else
            return x/(exp(x)-1)
        end
    end
    


    h=1.0/convert(Float64,n)
    geom=FVMGraph(collect(0:h:1))
    
    parameters=ILiqParameters()
    
    function flux!(this::ILiqParameters,f,uk,ul)
        f[iphi]=this.eps*(uk[iphi]-ul[iphi])
        muk=log(1-uk[ic])
        mul=log(1-ul[ic])

        #muk=0.0
        #mul=0.0
        f[ic]=bernoulli(this.a*(ul[iphi]-uk[iphi])+(muk-mul))*uk[ic]-bernoulli(this.a*(uk[iphi]-ul[iphi])+(mul-muk))*ul[ic]
    end 

    function storage!(this::FVMParameters, f,u)
        f[iphi]=0
        f[ic]=u[ic]
    end
    
    function reaction!(this::FVMParameters, f,u)
        f[iphi]=0.5-u[ic]
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
        inival[ic,inode]=0.5
    end
    parameters.eps=1.0e-3
    parameters.a=5
    control=FVMNewtonControl()
    control.verbose=true
    tstep=1.0e-2
    times=collect(0.0:tstep:1.0)
    for it=2:length(times)
        U=solve(sys,inival,control=control,tstep=tstep)
        for i=1:size(inival,2)
            inival[iphi,i]=U[iphi,i]
            inival[ic,i]=U[ic,i]
        end
        @printf("time=%g\n",times[it])
        if pyplot
            PyPlot.clf()
            plot(geom.Nodes[1,:],U[iphi,:])
            plot(geom.Nodes[1,:],U[ic,:])
        end
    end
    if pyplot && !isinteractive()
        pause(1.0e-10)
        waitforbuttonpress()
    end
end


if !isinteractive()
    run_iliq(args.n,pyplot=pyplot)
end




