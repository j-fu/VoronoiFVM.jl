#=

# 204: 2D Convection in Hagen-Poiseuille flow
([source code](SOURCE_URL))

Solve the equation

```math
\partial_t u -\nabla ( D \nabla u - v u) = 0
```
in $\Omega=(0,L)\times (0,H)$ with dirichlet boundary conditions
at $x=0$ and outflow boundary condition at $x=L$.
=#

module Example204_HagenPoiseuille
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize

function main(; nref = 0, Plotter = nothing, D = 0.01, v = 1.0, tend = 100, cin = 1.0)
    H = 1.0
    L = 5.0
    grid = simplexgrid(range(0, L; length = 20 * 2^nref),
                       range(0, H; length = 5 * 2^nref))

    function fhp(x, y)
        yh = y / H
        return v * 4 * yh * (1.0 - yh), 0
    end

    evelo = edgevelocities(grid, fhp)
    bfvelo = bfacevelocities(grid, fhp)

    function flux!(f, u, edge)
        vd = evelo[edge.index] / D
        bp = fbernoulli(vd)
        bm = fbernoulli(-vd)
        f[1] = D * (bp * u[1] - bm * u[2])
    end

    function outflow!(f, u, node)
        if node.region == 2
            f[1] = bfvelo[node.ibnode, node.ibface] * u[1]
        end
    end

    ispec = 1
    physics = VoronoiFVM.Physics(; flux = flux!, breaction = outflow!)
    sys = VoronoiFVM.DenseSystem(grid, physics)
    enable_species!(sys, ispec, [1])

    boundary_dirichlet!(sys, ispec, 4, cin)

    ## Transient solution of the problem
    control = VoronoiFVM.NewtonControl()
    control.Δt = 0.01 * 2.0^(-nref)
    control.Δt_min = 0.01 * 2.0^(-nref)
    control.Δt_max = 0.1 * tend
    control.force_first_step = true
    tsol = solve(sys; inival = 0, times = [0, tend], control = control)

    vis = GridVisualizer(; Plotter = Plotter)
    for i = 1:length(tsol.t)
        scalarplot!(vis[1, 1], grid, tsol[1, :, i]; flimits = (0, cin + 1.0e-5),
                    title = @sprintf("time=%3f", tsol.t[i]), show = true)
    end
    tsol
end

function test()
    tsol = main()
    all(tsol[end] .≈ 1)
end

end
