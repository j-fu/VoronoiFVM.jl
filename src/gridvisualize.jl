# Extend GridVisualize methods
import GridVisualize

"""
    $(TYPEDSIGNATURES)

Plot grid behind system
"""
GridVisualize.gridplot(sys::AbstractSystem; kwargs...) = GridVisualize.gridplot(sys.grid; kwargs...)

"""
    $(TYPEDSIGNATURES)

Plot grid behind system
"""
GridVisualize.gridplot!(vis, sys::AbstractSystem; kwargs...) = GridVisualize.gridplot!(vis, sys.grid; kwargs...)

"""
    $(TYPEDSIGNATURES)

Plot one species from solution
"""
function GridVisualize.scalarplot(sys::AbstractSystem, sol::AbstractMatrix; species = 1, scale=1.0, kwargs...)
    GridVisualize.scalarplot(sys.grid, sol[species, :]*scale; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Plot one species from solution
"""
function GridVisualize.scalarplot!(vis, sys::AbstractSystem, sol::AbstractMatrix; species = 1, scale=1.0, kwargs...)
    GridVisualize.scalarplot!(vis, sys.grid, sol[species, :]*scale; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Plot one species from transient solution
"""
function GridVisualize.scalarplot(sys::AbstractSystem, sol::AbstractTransientSolution; species = 1, scale=1.0, tscale=:identity, kwargs...)
    @assert dim_space(grid) == 1
    vis = GridVisualizer(kwargs...)
    if !isnothing(vis.Plotter)
        scalarplot!(vis[1, 1], sys, sol; species = species,scale=scale, tscale=tscale, kwargs...)
        reveal(vis)
    end
end

"""
    $(TYPEDSIGNATURES)

Plot one species from transient solution
"""
function GridVisualize.scalarplot!(vis, sys::AbstractSystem, sol::AbstractTransientSolution; species = 1,scale=1.0, tscale = :identity, tlabel = "t",
                                   kwargs...)
    if !isnothing(vis)
        grid = sys.grid
        @assert dim_space(grid) == 1
        f = vec(sol[species, :, :])*scale
        coord = grid[Coordinates]
        X = view(coord, 1, :)
        T = sol.t
        if tscale == :log
            T = log10.(sol.t)
            if sol.t[1] ≈ 0.0
                T[1] = log10(sol.t[2] / 2)
            end
            if tlabel == "t"
                tlabel = "log10(t)"
            end
        end
        (xmin, xmax) = extrema(X)
        (tmin, tmax) = extrema(T)
        xtaspect = 0.7 * (tmax - tmin) / (xmax - xmin)
        GridVisualize.scalarplot!(vis, simplexgrid(X, T), f; aspect = 1 / xtaspect, ylabel = tlabel, kwargs...)
    end
end

"""
    plothistory(sys, tsol; 
                plots=[:timesteps,:newtonsteps,:updates], 
                size=(700,400), 
                Plotter=GridVisualize.default_plotter(),
                kwargs...)

Plot solution history stored in `tsol`. The `plots` argument allows to choose which plots are shown.
"""
function plothistory(sys, tsol; plots=[:timesteps,:newtonsteps,:updates], size=(700,400), Plotter=GridVisualize.default_plotter(), kwargs...)
    hist=history(tsol)
    t=hist.times
    if length(t)==0
        error("Empty history. Did you pass `log=true` to the `solve()` method ?")
    end
    nplots=length(plots)
    vis=GridVisualize.GridVisualizer(;layout=(nplots,1), size, Plotter, kwargs...)

    for iplot in eachindex(plots)
        if plots[iplot]==:timesteps
            GridVisualize.scalarplot!(vis[iplot,1],t[1:(end - 1)], t[2:end] - t[1:(end - 1)];
                                      title = "Time step sizes",  xlabel = "t/s", ylabel = "Δt/s", color=:red, kwargs...)
        elseif plots[iplot]==:newtonsteps
            newtons=map(h->length(h.updatenorm),hist.histories)
            GridVisualize.scalarplot!(vis[iplot,1],t, newtons,xlabel = "t/s", ylabel = "n",
                        title="Newton iterations", color=:red, limits=(0,maximum(newtons)+1), kwargs...)
        elseif plots[iplot]==:updates
            GridVisualize.scalarplot!(vis[iplot,1],t, hist.updates,xlabel = "t/s", ylabel = "du",
                        title="Time step updates", color=:red, kwargs...)
        end
    end
    GridVisualize.reveal(vis)
end
