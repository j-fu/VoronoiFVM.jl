module Laplace

# We wrap this example and all later ones
# into a module structure. This allows to load
# all of them at once into the REPL without name
# clashes. We shouldn't forget the corresponding end
# statement.
using VoronoiFVM



# Flux function which describes the flux
# between neigboring control volumes
function g(f,u,edge,data)
    uk=viewK(edge,u)
    ul=viewL(edge,u)
    f[1]=uk[1]-ul[1]
end


function main()
    
    nspecies=1
    ispec=1

    # Create a one dimensional discretization project
    X=collect(0:0.2:1)
    grid=VoronoiFVM.Grid(X)


    # Create a physics structure
    physics=VoronoiFVM.Physics(num_species=nspecies,flux=g)

    # Create a finite volume system
    sys=VoronoiFVM.DenseSystem(grid,physics)

    # Enable species 1 in region 1
    enable_species!(sys,ispec,[1])

    # Set boundary conditions
    sys.boundary_factors[ispec,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[ispec,2]=VoronoiFVM.Dirichlet
    sys.boundary_values[ispec,1]=0
    sys.boundary_values[ispec,2]=1
    
    # Create & initialize array for solution and initial value
    inival=unknowns(sys)
    solution=unknowns(sys)
    inival.=0
    solution.=0

    # Solve stationary problem
    solve!(solution,inival,sys)

    # Return test value
    return solution[3]
end



end
