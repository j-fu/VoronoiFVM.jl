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
