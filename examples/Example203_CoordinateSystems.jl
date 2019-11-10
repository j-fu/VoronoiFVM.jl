#= 

# 203: Various coordinate systems
=#

module Example203_CoordinateSystems

using VoronoiFVM
using LinearAlgebra
using Plots


function main(;nref=0,r1=0.0, r2=5.0, dim=2,doplot=false)
    h=0.1*2.0^(-nref)
    R=collect(r1:h:r2)
    Z=collect(0:h:2)
    grid=VoronoiFVM.Grid(R)
    if dim==2
        grid=VoronoiFVM.Grid(R,Z)
    end
    cylindrical!(grid)
    
    function symlapcyl(r)
        if r1≈0.0
            return 0.25*(r2*r2-r*r);
        else
            return (log(r)-log(r2))/(log(r1)-log(r2));
        end
    end

    function flux!(f,u,edge,data)
        uk=viewK(edge,u)  
        ul=viewL(edge,u)
        f[1]=uk[1]-ul[1]
    end
    
    function source!(f,node,data)
        if r1≈0.0
            f[1]=1.0
        else
            f[1]=0.0
        end
    end
    
    # Create a physics structure
    physics=VoronoiFVM.Physics(num_species=1,flux=flux!,source=source!)
    sys=VoronoiFVM.DenseSystem(grid,physics)
    ispec=1
    enable_species!(sys,ispec,[1])
    ileft=1
    if dim==2
        ileft=4
    end
    if !(r1≈0.0)
        sys.boundary_factors[ispec,ileft]=VoronoiFVM.Dirichlet
        sys.boundary_values[ispec,ileft]=1
    end
    sys.boundary_factors[ispec,2]=VoronoiFVM.Dirichlet
    sys.boundary_values[ispec,2]=0
    
    inival=unknowns(sys)
    solution=unknowns(sys)
    inival.=0
    solution.=0

    # Solve stationary problem
    solve!(solution,inival,sys)
    
    if doplot
        p=contourf(R,Z,transpose(reshape(values(solution),length(R),length(Z))),colorbar=:right)
        gui(p)
    end
    
    exact=symlapcyl.(grid.coord[1,:])
    err=norm(solution[1,:]-exact,Inf)
end

#
# Called by unit test
#
function test()
    main() <1.0e-14
end

end
