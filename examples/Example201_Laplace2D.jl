#=

# 201: 2D Laplace equation 
([source code](SOURCE_URL))


=#

module Example201_Laplace2D

using VoronoiFVM,ExtendableGrids


## Flux function which describes the flux
## between neigboring control volumes
function g!(f,u0,edge)
    u=unknowns(edge,u0)
    f[1]=u[1,1]-u[1,2]
end


function main(;Plotter=nothing)
    nspecies=1 
    ispec=1    
    X=collect(0:0.2:1)
    grid=VoronoiFVM.Grid(X,X)
    physics=VoronoiFVM.Physics(num_species=nspecies,flux=g!)
    sys=VoronoiFVM.System(grid,physics)
    enable_species!(sys,ispec,[1])
    boundary_dirichlet!(sys,ispec,1,0.0)
    boundary_dirichlet!(sys,ispec,3,1.0)
    inival=unknowns(sys,inival=0)
    solution=unknowns(sys)
    solve!(solution,inival,sys)
    gridplot(grid,solution[1,:],Plotter=Plotter)
    return solution[7]
end

## Called by unit test

function test()
    main() ≈ 0.2
end

end

