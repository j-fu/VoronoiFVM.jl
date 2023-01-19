#=

# 160: Unipolar degenerate drift-diffusion
([source code](SOURCE_URL))

See: C. Cancès, C. Chainais-Hillairet, J. Fuhrmann, and B. Gaudeul, "A numerical-analysis-focused comparison of several finite volume schemes for a unipolar degenerate drift-diffusion model" IMA Journal of Numerical Analysis, vol. 41, no. 1, pp. 271–314, 2021.  

Available from [https://doi.org/10.1093/imanum/draa002](https://doi.org/10.1093/imanum/draa002),
the preprint is on [arxiv1907.11126](https://arxiv.org/abs/1907.11126).

The problem consists of a Poisson equation for the electrostatic potential $\phi$:

```math
-\nabla \varepsilon \nabla \phi = z(2c-1)
```
and a degenerate drift-diffusion equation of the transport of a charged species $c$:

```math
\partial_t u  - \nabla\cdot \left(\nabla c  + c \nabla (\phi - \log (1-c) )\right)
```

In particular, the paper, among others, investigates the "sedan" flux discretization which is able to handle the degeneracy coming from the $\log (1-c)$ term. The earliest reference to this scheme we found in the source code of the [SEDAN III](http://www-tcad.stanford.edu/tcad/programs/sedan3.html) semiconductor device simulator.
=#

module Example160_UnipolarDriftDiffusion1D

using Printf

using VoronoiFVM
using ExtendableGrids
using GridVisualize

mutable struct Data
    eps::Float64
    z::Float64
    ic::Int32
    iphi::Int32
    Data() = new()
end

function classflux!(f, u, edge, data)
    ic = data.ic
    iphi = data.iphi
    f[iphi] = data.eps * (u[iphi, 1] - u[iphi, 2])
    bp, bm = fbernoulli_pm(u[iphi, 1] - u[iphi, 2])
    f[ic] = bm * u[ic, 1] - bp * u[ic, 2]
end

function storage!(f, u, node, data)
    ic = data.ic
    iphi = data.iphi
    f[iphi] = 0
    f[ic] = u[ic]
end

function reaction!(f, u, node, data)
    ic = data.ic
    iphi = data.iphi
    f[iphi] = data.z * (1 - 2 * u[ic])
    f[ic] = 0
end

function sedanflux!(f, u, edge, data)
    ic = data.ic
    iphi = data.iphi
    f[iphi] = data.eps * (u[iphi, 1] - u[iphi, 2])
    mu1 = -log(1 - u[ic, 1])
    mu2 = -log(1 - u[ic, 2])
    bp, bm = fbernoulli_pm(data.z * 2 * (u[iphi, 1] - u[iphi, 2]) + (mu1 - mu2))
    f[ic] = bm * u[ic, 1] - bp * u[ic, 2]
end

function main(; n = 20, Plotter = nothing, dlcap = false, verbose = false,
              unknown_storage = :sparse, DiffEq = nothing)
    h = 1.0 / convert(Float64, n)
    grid = VoronoiFVM.Grid(collect(0:h:1))

    data = Data()
    data.eps = 1.0e-3
    data.z = -1
    data.iphi = 1
    data.ic = 2

    ic = data.ic
    iphi = data.iphi

    physics = VoronoiFVM.Physics(; data = data,
                                 flux = sedanflux!,
                                 reaction = reaction!,
                                 storage = storage!)
    sys = VoronoiFVM.System(grid, physics; unknown_storage = unknown_storage)

    enable_species!(sys, 1, [1])
    enable_species!(sys, 2, [1])

    boundary_dirichlet!(sys, iphi, 1, 5.0)
    boundary_dirichlet!(sys, iphi, 2, 0.0)
    boundary_dirichlet!(sys, ic, 2, 0.5)

    inival = unknowns(sys)
    @views inival[iphi, :] .= 2
    @views inival[ic, :] .= 0.5

    if !dlcap
        ## Create solver control info for constant time step size
        tstep = 1.0e-3
        control = VoronoiFVM.NewtonControl()
        control.verbose = verbose
        control.Δt_min = tstep
        control.Δt = tstep
        control.Δt_grow = 1.2
        control.Δt_max = 0.1
        control.Δu_opt = 100
        control.damp_initial = 0.5
        if isnothing(DiffEq)
            tsol = solve(inival, sys, [0.0, 10]; control = control)
        else # does not work yet...
            tsol = solve(DiffEq, inival, sys, [0.0, 10];
                         initializealg = DiffEq.NoInit(),
                         dt = tstep)
        end
        vis = GridVisualizer(; Plotter = Plotter, layout = (1, 1), fast = true)
        for log10t = -4:0.01:0
            time = 10^(log10t)
            sol = tsol(time)
            scalarplot!(vis[1, 1], grid, sol[iphi, :]; label = "ϕ",
                        title = @sprintf("time=%.3g", time), flimits = (0, 5),
                        color = :green)
            scalarplot!(vis[1, 1], grid, sol[ic, :]; label = "c", flimits = (0, 5),
                        clear = false, color = :red)
            reveal(vis)
        end
        return sum(tsol[end])

    else  # Calculate double layer capacitance
        U = unknowns(sys)
        control = VoronoiFVM.NewtonControl()
        control.damp_initial = 1.0e-5
        delta = 1.0e-4
        @views inival[iphi, :] .= 0
        @views inival[ic, :] .= 0.5
        sys.boundary_values[iphi, 1] = 0

        dphi = 1.0e-1
        phimax = 5
        delta = 1.0e-4
        vplus = zeros(0)
        cdlplus = zeros(0)
        vminus = zeros(0)
        cdlminus = zeros(0)
        cdl = 0.0
        vis = GridVisualizer(; Plotter = Plotter, layout = (2, 1), fast = true)
        for dir in [1, -1]
            phi = 0.0
            while phi < phimax
                sys.boundary_values[iphi, 1] = dir * phi
                solve!(U, inival, sys; control = control)
                inival .= U
                Q = integrate(sys, physics.reaction, U)
                sys.boundary_values[iphi, 1] = dir * phi + delta
                solve!(U, inival, sys; control = control)
                inival .= U

                scalarplot!(vis[1, 1], grid, U[iphi, :]; label = "ϕ",
                            title = @sprintf("Δϕ=%.3g", phi), flimits = (-5, 5),
                            clear = true, color = :green)
                scalarplot!(vis[1, 1], grid, U[ic, :]; label = "c", flimits = (0, 5),
                            clear = false, color = :red)

                Qdelta = integrate(sys, physics.reaction, U)
                cdl = (Qdelta[iphi] - Q[iphi]) / delta
                if dir == 1
                    push!(vplus, dir * phi)
                    push!(cdlplus, cdl)
                else
                    push!(vminus, dir * phi)
                    push!(cdlminus, cdl)
                end

                scalarplot!(vis[2, 1], [0, 1.0e-1], [0, 0.05]; color = :white, clear = true)
                v = vcat(reverse(vminus), vplus)
                c = vcat(reverse(cdlminus), cdlplus)
                if length(v) >= 2
                    scalarplot!(vis[2, 1], v, c; color = :green, clear = false,
                                title = "C_dl")
                end

                phi += dphi
                reveal(vis)
            end
        end

        return cdl
    end
end

function test()
    isapprox(main(; unknown_storage = :sparse, dlcap = false), 18.721369939561963;
             rtol = 1.0e-5) &&
        isapprox(main(; unknown_storage = :sparse, dlcap = true), 0.010759276468375045;
                 rtol = 1.0e-5) &&
        isapprox(main(; unknown_storage = :dense, dlcap = false), 18.721369939561963;
                 rtol = 1.0e-5) &&
        isapprox(main(; unknown_storage = :dense, dlcap = true), 0.010759276468375045;
                 rtol = 1.0e-5)
end
end
