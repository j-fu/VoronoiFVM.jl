# Extend GridVisualize methods
import GridVisualize

GridVisualize.gridplot(sys::System; kwargs...)=GridVisualize.gridplot(sys.grid; kwargs...)

GridVisualize.gridplot!(vis,sys::System; kwargs...)=GridVisualize.gridplot!(vis,sys.grid; kwargs...)

function GridVisualize.scalarplot(sys::System, sol::AbstractMatrix; species=1, kwargs...)
    GridVisualize.scalarplot(sys.grid,sol[species,:]; kwargs...)
end

function GridVisualize.scalarplot!(vis,sys::System, sol::AbstractVector; species=1, kwargs...)
    GridVisualize.scalarplot!(vis,sys.grid,sol[species,:], kwargs...)
end
