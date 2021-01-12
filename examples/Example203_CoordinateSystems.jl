#=

# 203: Various coordinate systems 
([source code](SOURCE_URL))
=#

module Example203_CoordinateSystems

using VoronoiFVM
using LinearAlgebra
using ExtendableGrids
using GridVisualize

function main(;nref=0,r1=0.0, r2=5.0, dim=2,Plotter=nothing)
    h=0.1*2.0^(-nref)
    R=collect(r1:h:r2)
    Z=collect(0:h:2)
    grid=VoronoiFVM.Grid(R)
    if dim==2
        grid=VoronoiFVM.Grid(R,Z)
    end
    circular_symmetric!(grid)
    
    function symlapcyl(r)
        if r1≈0.0
            return 0.25*(r2*r2-r*r);
        else
            return (log(r)-log(r2))/(log(r1)-log(r2));
        end
    end

    function flux!(f,u0,edge)
        u=unknowns(edge,u0)
        f[1]=u[1,1]-u[1,2]
    end
    
    function source!(f,node)
        if r1≈0.0
            f[1]=1.0
        else
            f[1]=0.0
        end
    end
    
    # Create a physics structure
    physics=VoronoiFVM.Physics(num_species=1,flux=flux!,source=source!)
    sys=VoronoiFVM.System(grid,physics,unknown_storage=:dense)
    ispec=1
    enable_species!(sys,ispec,[1])
    ileft=1
    if dim==2
        ileft=4
    end
    if !(r1≈0.0)
        boundary_dirichlet!(sys,ispec,ileft,1.0)
    end
    boundary_dirichlet!(sys,ispec,2,0.0)
     
    inival=unknowns(sys)
    solution=unknowns(sys)
    inival.=0
    solution.=0

    # Solve stationary problem
    solve!(solution,inival,sys)

    scalarplot(grid,solution[1,:],Plotter=Plotter)

    exact=symlapcyl.(coordinates(grid)[1,:])
    err=norm(solution[1,:]-exact,Inf)
end

#
# Called by unit test
#
function test()
    main() <1.0e-14
end

end
