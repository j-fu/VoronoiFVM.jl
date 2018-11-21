module OneSpeciesNonlinearPoisson

using Printf
using TwoPointFluxFVM

if isinteractive()
    using PyPlot
end


"""
Structure containing  userdata information
"""
mutable struct Physics <:TwoPointFluxFVM.Physics
    TwoPointFluxFVM.@AddPhysicsBaseClassFields
    eps::Float64 
    Physics()=Physics(new())
end

"""
Reaction term
"""
function reaction!(this::Physics,f::AbstractArray,u::AbstractArray)
    f[1]=u[1]^2
end

"""
Flux term
"""
function flux!(this::Physics,f::AbstractArray,uk::AbstractArray,ul::AbstractArray)
    f[1]=this.eps*(uk[1]^2-ul[1]^2)
end 


"""
Source term
"""
function source!(this::Physics,f::AbstractArray,x::AbstractArray)
    f[1]=1.0e-4*x[1]
end 

"""
Constructor for userdata structure
"""
function Physics(this::Physics)
    TwoPointFluxFVM.PhysicsBase(this,1)
    this.eps=1
    this.flux=flux!
    this.reaction=reaction!
    this.source=source!
    return this
end

"""
Main function for user interaction from REPL and
for test. Default parameters need to generate correct
test value.
"""
function main(;n=10,pyplot=false,verbose=false)
    h=1.0/convert(Float64,n)
    geom=TwoPointFluxFVM.Graph(collect(0:h:1))

    physics=Physics()

    sys=TwoPointFluxFVM.System(geom,physics)
    sys.boundary_values[1,1]=1.0
    sys.boundary_values[1,2]=0.5
    
    sys.boundary_factors[1,1]=TwoPointFluxFVM.Dirichlet
    sys.boundary_factors[1,2]=TwoPointFluxFVM.Dirichlet
    
    inival=unknowns(sys)
    inival.=0.5

    physics.eps=1.0e-2

    control=TwoPointFluxFVM.NewtonControl()
    control.verbose=verbose
    tstep=1.0e-2
    times=collect(0.0:tstep:1.0)
    test_result=0
    for it=2:length(times)
        U=solve(sys,inival,control=control,tstep=tstep)
        test_result=U[5]
        for i in eachindex(U)
            inival[i]=U[i]
        end
        if verbose
            @printf("time=%g\n",times[it])
        end
        if pyplot
            PyPlot.clf()
            PyPlot.grid()
            plot(geom.node_coordinates[1,:],U)
            pause(1.0e-10)
        end
    end
    return test_result
end

end

