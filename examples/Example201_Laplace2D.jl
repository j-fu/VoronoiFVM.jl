#=

# 201: 2D Laplace equation 
([source code](SOURCE_URL))


=#

module Example201_Laplace2D

using VoronoiFVM,ExtendableGrids
using GridVisualize

## Flux function which describes the flux
## between neigboring control volumes
function g!(f,u,edge)
    f[1]=u[1,1]-u[1,2]
end


function main(;Plotter=nothing,n=5)
    nspecies=1 
    ispec=1    
    X=collect(0:1.0/n:1)
    grid=VoronoiFVM.Grid(X,X)


    physics=VoronoiFVM.Physics(flux=g!)
    sys=VoronoiFVM.System(grid,physics)
    enable_species!(sys,ispec,[1])
    boundary_dirichlet!(sys,ispec,1,0.0)
    boundary_dirichlet!(sys,ispec,3,1.0)
    inival=unknowns(sys,inival=0)
    solution=unknowns(sys)
    solve!(solution,inival,sys)
    nf=nodeflux(sys,solution)
    vis=GridVisualizer(Plotter=Plotter)
    scalarplot!(vis,grid,solution[1,:],clear=true,colormap=:summer)
    vectorplot!(vis,grid,nf[:,1,:],clear=false,spacing=0.1,vscale=0.5)
    reveal(vis)
    return solution[7]
end

## Called by unit test

function test()
    main() ≈ 0.2
end

end

