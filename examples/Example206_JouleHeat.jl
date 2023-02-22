# # 206: 2D Joule heating
# ([source code](SOURCE_URL))
#
## |------------------------------------------------------------------|
## | Joule Heat System						   |
## |	 - div (kappa  grad phi)   = 0				   |
## | dt (cT)  - div (lambda grad T)   = kappa (grad phi) * (grad phi) |
## |			     kappa = kappa0 exp(alpha (T-T0))	   |
## |------------------------------------------------------------------|


module Example206_JouleHeat

using Printf
using VoronoiFVM
using ExtendableGrids
using ExtendableSparse
using GridVisualize
using LinearAlgebra
using SimplexGridFactory
using Triangulate

function main(; nref = 0, Plotter = nothing, verbose = "and", unknown_storage = :sparse,     ythin=0.25)


    b=SimplexGridBuilder(Generator=Triangulate)
    p00=point!(b,0,0)
    p30=point!(b,3,0)
    p32=point!(b,3,1)
    p21=point!(b,2,ythin)
    p11=point!(b,1,ythin)
    p02=point!(b,0,1)

    facetregion!(b,4)
    facet!(b,p00,p30)
    facetregion!(b,2)
    facet!(b,p30,p32)
    facetregion!(b,3)
    facet!(b,p32,p21)
    facet!(b,p21,p11)
    facet!(b,p11,p02)
    facetregion!(b,1)
    facet!(b,p02,p00)
    
    grid=simplexgrid(b;maxvolume=0.01*4.0^(-nref))
    
    iϕ::Int=1
    iT::Int=2
    κ0::Float64=1;
    α::Float64=1;
    T0::Float64=0.5;
    λ::Float64=1;
    c::Float64=1;
    
    function storage!(y,u,node)
        y[iT]=c*u[iT]
    end
    
    κ(T)=κ0*exp(α*(T-T0))
    
    function flux!(y,u,edge)
        y[iϕ]= κ(y[iT])*(u[iϕ,1]-u[iϕ,2])
        y[iT]= λ*(u[iT,1]-u[iT,2])
    end
    

    function jouleheat!(y,u,edge)
        y[iT]= -κ(y[iT])*(u[iϕ,1]-u[iϕ,2])*(u[iϕ,1]-u[iϕ,2])
    end
    
    function bcondition!(y,u,node)
        boundary_dirichlet!(y,u,node,species=iϕ,region=1, value=-10)
        boundary_dirichlet!(y,u,node,species=iϕ,region=2, value=10)

        boundary_robin!(y,u,node,species=iT,region=1, value=T0, factor=0.5)
        boundary_robin!(y,u,node,species=iT,region=2, value=T0, factor=0.5)
        boundary_robin!(y,u,node,species=iT,region=3, value=T0, factor=0.5)
        boundary_robin!(y,u,node,species=iT,region=4, value=T0, factor=0.5)
        
    end
    
    sys=VoronoiFVM.System(grid; bcondition=bcondition!, flux=flux!, edgereaction=jouleheat!,storage=storage!,species=[iϕ,iT])

    sol=solve(sys;verbose)

    vis=GridVisualizer(;Plotter,layout=(2,1))
    scalarplot!(vis[1,1],grid,sol[iϕ,:],title="ϕ",colormap=:bwr)
    scalarplot!(vis[2,1],grid,sol[iT,:],title="T", colormap=:hot)
    reveal(vis)
    norm(sol,Inf)
end

function test()
    testval=24.639120035942938
    main() ≈ testval
end
end
