module OneSpeciesNonlinearPoisson

# This gives us he @printf macro (c-like output)
using Printf

# That's the thing we want to do
using VoronoiFVM

# Allow plotting
if isinteractive()
    using PyPlot
end




# Main function for user interaction from REPL and
# for testimg. Default physics need to generate correct
# test value.

function main(;n=10,pyplot=false,verbose=false, dense=false)



    # Create a one dimensional discretization project
    h=1.0/convert(Float64,n)
    grid=VoronoiFVM.Grid(collect(0:h:1))

    eps=1.0e-2
    

    # Flux function which describes the flux
    # between neigboring control volumes
    flux=function(f,u,edge,data)
        f[1]=eps*(u[1]^2-u[2]^2)
    end


    # Source term
    source=function(f,node,data)
        f[1]=1.0e-4*node.coord[1]
    end

    # Storage term (under the time derivative)
    storage=function(f,u,node,data)
        f[1]=u[1]
    end
    
    # Reation term
    reaction=function(f,u,node,data)
        f[1]=u[1]^2
    end
    # Create a physics structure
    physics=VoronoiFVM.Physics(flux=flux,source=source,storage=storage,reaction=reaction)


    # Create a finite volume system - either
    # in the dense or  the sparse version.
    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics)
    end

    # Add species 1 to region 1
    enable_species!(sys,1,[1])

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
        solve!(U,inival,sys, control=control,tstep=tstep)
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

