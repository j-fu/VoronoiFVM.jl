# We wrap this example and all later ones
# into a module structure. This allows to load
# all of them at once into the REPL without name
# clashes. We shouldn't forget the corresponding end
# statement.
module OneSpeciesNonlinearPoisson

# This gives us he @printf macro (c-like output)
using Printf

# That's the thing we want to do
using VoronoiFVM

# Allow plotting
if isinteractive()
    using PyPlot
end



#
# Structure containing  userdata information
#
# We choose a mutable struct which allows to overwrite
# fields later.
mutable struct Physics  <: VoronoiFVM.Physics
    flux::Function      # flux function, mandatory
    source::Function    # source function, optional
    reaction::Function  # reaction term, optional
    storage::Function   # storage term, mandatory
    eps::Real           # Example for "user data" passed to the callback

    Physics()=new() # Provide inner constructor resulting in uninitialized struct
end



# Main function for user interaction from REPL and
# for testimg. Default physics need to generate correct
# test value.

function main(;n=10,pyplot=false,verbose=false, dense=false)



    # Create a one dimensional discretization project
    h=1.0/convert(Float64,n)
    grid=VoronoiFVM.Grid(collect(0:h:1))
    
    # Create a physics structure
    physics=Physics()
    physics.eps=1.0e-2

    # Flux function which describes the flux
    # between neigboring control volumes
    physics.flux=function(physics,edge,f,u)
        nspecies=1
        uk=VoronoiFVM.UK(u,nspecies)
        ul=VoronoiFVM.UL(u,nspecies)
        f[1]=physics.eps*(uk[1]^2-ul[1]^2)
    end 


    # Source term
    physics.source=function(physics,node,f)
        f[1]=1.0e-4*node.coord[1]
    end 

    # Storage term (under the time derivative)
    physics.storage=function(physics,node, f,u)
        f[1]=u[1]
    end

    # Reation term
    physics.reaction=function(physics,node,f,u)
        f[1]=u[1]^2
    end


    # Create a finite volume system - either
    # in the dense or  the sparse version.
    # Need to provide the overall number of species here
    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics,1)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics,1)
    end

    # Add species 1 to region 1
    add_species(sys,1,[1])

    # Set boundary conditions
    sys.boundary_values[1,1]=1.0
    sys.boundary_values[1,2]=0.5
    sys.boundary_factors[1,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[1,2]=VoronoiFVM.Dirichlet

    # Create a solution array
    inival=unknowns(sys)
    U=unknowns(sys)

    # Broadcast the initial value
    inival.=0.5

    # Create solver control info
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose

    # time stepping
    tstep=1.0e-2
    times=collect(0.0:tstep:1.0)
    test_result=0
    for it=2:length(times)
        # Solve for new timestep with old timestep
        # solution in inival
        solve(sys,inival,U, control=control,tstep=tstep)
        test_result=U[5]

        # Update inival
        inival.=U

        if verbose
            @printf("time=%g\n",times[it])
        end

        # Plot data
        if pyplot
            PyPlot.clf()
            PyPlot.grid()
            plot(grid.coord[1,:],U[1,:])
            pause(1.0e-10)
        end
    end
    return test_result
end


end # Yes, this is *that* end (of the module...)

