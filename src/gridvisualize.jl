# Extend GridVisualize methods
import GridVisualize

GridVisualize.gridplot(sys::AbstractSystem; kwargs...)=GridVisualize.gridplot(sys.grid; kwargs...)

GridVisualize.gridplot!(vis,sys::AbstractSystem; kwargs...)=GridVisualize.gridplot!(vis,sys.grid; kwargs...)

function GridVisualize.scalarplot(sys::AbstractSystem, sol::AbstractMatrix; species=1, kwargs...)
    GridVisualize.scalarplot(sys.grid,sol[species,:]; kwargs...)
end

function GridVisualize.scalarplot!(vis,sys::AbstractSystem, sol::AbstractMatrix; species=1, kwargs...)
    GridVisualize.scalarplot!(vis,sys.grid,sol[species,:]; kwargs...)
end


function GridVisualize.scalarplot(sys::AbstractSystem, sol::TransientSolution; species=1, kwargs...)
    @assert dim_space(grid)==1
    vis=GridVisualizer(kwargs...)
    if !isnothing(vis.Plotter)
        scalarplot!(vis[1,1],sys,sol; species=species, kwargs...)
        reveal(vis)
    end
end



function GridVisualize.scalarplot!(vis,sys::AbstractSystem, sol::TransientSolution; species=1, kwargs...)
    if !isnothing(vis)
        grid=sys.grid
        @assert dim_space(grid)==1
        f=vec(sol[species,:,:])
        coord=grid[Coordinates]
        GridVisualize.scalarplot!(vis,simplexgrid(coord[1,:],sol.t),f; kwargs...)
    end
end
