#=

# 301: 3D Laplace equation 
([source code](SOURCE_URL))


=#

module Example301_Laplace3D

using VoronoiFVM,ExtendableGrids
using GridVisualize


## Flux function which describes the flux
## between neigboring control volumes
function g!(f,u0,edge)
    u=unknowns(edge,u0)
    f[1]=u[1,1]-u[1,2]
end

function s(f,node)
    n=view(node.coord,:,node.index)
    f[1]=n[1]*sin(5.0*n[2])*exp(n[3])
end


function main(;Plotter=nothing,n=5)
    nspecies=1 
    ispec=1    
    X=collect(0:1/n:1)
    grid=VoronoiFVM.Grid(X,X,X)
    physics=VoronoiFVM.Physics(flux=g!,source=s)
    sys=VoronoiFVM.System(grid,physics)
    enable_species!(sys,ispec,[1])
    boundary_dirichlet!(sys,ispec,5,0.0)
    boundary_dirichlet!(sys,ispec,6,0.0)
    inival=unknowns(sys,inival=0)
    solution=unknowns(sys)
    solve!(solution,inival,sys)
    scalarplot(grid,solution[1,:],Plotter=Plotter,zplane=0.5, flevel=0.5)
    return solution[43]
end

## Called by unit test

function test()
    main() ≈ 0.012234524449380824 
end

end

