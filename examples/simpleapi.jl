module ExampleSimpleAPI

using VoronoiFVM
using GridVisualize
using ExtendableGrids



function xmain(;Plotter=nothing,n=50)
    X=0:1.0/n:1
    grid=simplexgrid(X,X)

    flux(u,edge)= u[1,1]-u[1,2]
    reaction(u,node)=  u[1,1]^3-10

    bc=[
        DirichletBC(species=1,boundary=1,value=0),
        DirichletBC(species=1,boundary=3,value=1),
    ]
    
    sys=system(grid; num_species=1,flux=flux,reaction=reaction,bc=bc)
    solution=solve(sys; verbose=true, damp_initial=0.01)
    scalarplot(grid,solution[1,:]; Plotter=Plotter, clear=true,colormap=:summer,show=true)
    return sum(solution)
end


function oldmain(;Plotter=nothing,n=50)
    nspecies=1 
    ispec=1    
    X=collect(0:1.0/n:1)
    grid=VoronoiFVM.Grid(X,X)

    function flux(y,u,edge)
        y[1]=u[1,1]-u[1,2]
    end
    function reaction(y,u,node)
        y[1]=u[1,1]^3-10.0
    end
    
    physics=VoronoiFVM.Physics(flux=flux,reaction=reaction)
    sys=VoronoiFVM.System(grid,physics)
    enable_species!(sys,ispec,[1])
    boundary_dirichlet!(sys,ispec,1,0.0)
    boundary_dirichlet!(sys,ispec,3,1.0)
    inival=unknowns(sys,inival=0)
    solution=unknowns(sys)
    control=VoronoiFVM.NewtonControl()

    solve!(solution,inival,sys,control=control)
    scalarplot(grid,solution[1,:]; Plotter=Plotter, clear=true,colormap=:summer,show=true)
    return sum(solution)
end


## Called by unit test

function test()
    main() ≈ 0.2
end

end

