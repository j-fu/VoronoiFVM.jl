module Example000_SimpleAPI

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

Test callbacks on system creation (e.g. against allocations)

Simplifying changes: kwargs in solver, system creation
Improved interface: bc now in breaction callback.



struct BoundaryCondition
    ispec::Int32
    ireg::Int32
    fac::Float64
    val::Float64
end


neumann(;species=1, region=1,value=0)=BoundaryCondition(species,region,0,value)
robin(;species=1, region=1,factor=0,value=0)=BoundaryCondition(species,region,factor,value)
dirichlet(;species=1, region=1,value=0)=BoundaryCondition(species,region,Dirichlet,value)

function boundary_conditions!(y,u,bnode,bc)
    for bc in bconds
        if bc.fac==Dirichlet
            boundary_dirichlet!(y,u,bnode,bc.ispec,bc.ireg,bc.val)
        else
            boundary_robin!(y,u,bnode,bc.ispec,bc.ireg,bc.fac,bc.val)
        end
    end
end

function boundary_conditions(bconds)
    (y,u,bnode,args...) -> boundary_conditions!(y,u,bnode,bc)
end

    bc=boundary_conditions([dirichlet(region=1,value=0),
                            dirichlet(region=3,value=1)])

KISS 


"""
function main(;Plotter=nothing,n=50)
    X=0:1.0/n:1
    grid=simplexgrid(X,X)
    ispec=1

    
    flux!(y,u,edge)= y[1]= u[1,1]-u[1,2]

    reaction!(y,u,node)=  y[1]= u[1,1]^3-10
    
    
    function bc!(args...)
        boundary_dirichlet!(args...,region=1,value=0)
        boundary_dirichlet!(args...,region=3,value=1)
    end
    
    sys=VoronoiFVM.System(grid; species=1,flux=flux!,reaction=reaction!,breaction=bc!)
    solution=solve(sys; verbose=false, blob=true)
    scalarplot(grid,solution[1,:]; Plotter=Plotter, clear=true,colormap=:summer,show=true)
    return sum(solution)
end



function tmain(;Plotter=nothing,n=50)
    X=0:1.0/n:1
    grid=simplexgrid(X,X)
    ispec=1

    
    flux!(y,u,edge)= y[1]= u[1,1]-u[1,2]

    reaction!(y,u,node)=  y[1]= u[1,1]^3-10

    storage!(y,u,node)= y.=u
    
    function bc!(y,u,bnode)
        boundary_dirichlet!(y,u,bnode,region=1,value=0)
        boundary_dirichlet!(y,u,bnode,region=3,value=ramp(bnode.time; dt=(0,1.0e-2), du=(0,1) ))
    end
    
    sys=VoronoiFVM.System(grid; inival=0, species=1,storage=storage!,flux=flux!,reaction=reaction!,breaction=bc!)
    solution=solve(sys;times=[0,1], verbose=true)
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
    main() ≈ 2935.324529621136
end



end

