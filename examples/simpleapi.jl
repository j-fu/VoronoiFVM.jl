module ExampleSimpleAPI

using VoronoiFVM
using GridVisualize
using ExtendableGrids

"""

- Keep the mutating functions, just document them better and  provide tools for standard cases
- Develop patterns for wrapper functions
- This solves parameter handling for BC and the case of nonatonomous systems! 
   t=time(node)
   node.time
   or pass t as a parameter: todo: square this with the impedance stuff.
    
- Find pattern for parameter, species, quantity naming
- For natural BC no problem: just use breaction
- For Dirichlet:  Extra function type ? No, just give a pattern how to do this.
   function breaction!(y,u,bnode,data)
        boundary_dirichlet!(y,bnode,ispec,ireg,val)
        boundary_neumann!(y,bnode,ispec,ireg,val)
        boundary_robin!(y,bnode,ispec,ireg,fac,val)
   end
 
   This will allow at once to handle impedance for drift-diffusion and space + time + parameter dependent bc!

   General impedance approach just uses a parameter for excitation, and we will 
   get general impedance handling from this.
- Provide tools for building elementary flux functions: diffusion, upwind, scharfetter_gummel
- Work with species

Reply to rveltz:
Constraints: 
- callbacks able to work without allocations and with forwarddiff
- species numbers etc. at runtime => difficult to work with static arrays etc. in API
- => Mutating callbacks

Simplifying changes: kwargs in solver, system creation
Improved interface: bc now in breaction callback.
"""

function xmain(;Plotter=nothing,n=50)
    X=0:1.0/n:1
    grid=simplexgrid(X,X)
    ispec=1

    
    flux!(y,u,edge)= y[1]= u[1,1]-u[1,2]
    reaction!(y,u,node)=  y[1]= u[1,1]^3-10
    function bc!(args...)
        boundary_dirichlet!(args...,region=1,value=0)
        boundary_dirichlet!(args...,region=3,value=1)
    end
    
    sys=system(grid; num_species=1,flux=flux!,reaction=reaction!,breaction=bc!)
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

#=
=#

end
b
