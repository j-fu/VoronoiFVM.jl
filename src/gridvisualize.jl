# Extend GridVisualize methods


GridVisualize.gridplot(sys::System; kwargs...)=GridVisualize.gridplot(sys.grid; kwargs...)
GridVisualize.gridplot!(vis,sys::System; kwargs...)=GridVisualize.gridplot!(vis,sys.grid; kwargs...)


function GridVisualize.scalarplot(sys::System, sol::AbstractMatrix; kwargs...)
    nspec=size(sol,1)
    if nspec==1
        GridVisualize.scalarplot(sys.grid,sol[1,:]; kwargs...)
    end
end

function GridVisualize.scalarplot!(vis,sys::System, sol::AbstractMatrix; labels=nothing, kwargs...)
    nspec=size(sol,1)
    if nspec==1
        GridVisualize.scalarplot!(vis,sys.grid,sol[1,:]; kwargs...)
    elseif dim_space(sys.grid)==1
        cmap=GridVisualize.bregion_cmap(nspec)
        if isnothing(labels)
            labels=[@sprintf("u%d",ispec) for ispec=1:nspec]
        end
        clear=true
        for ispec=1:nspec
            GridVisualize.scalarplot!(vis,sys.grid,sol[ispec,:]; clear=clear,color=cmap[ispec], label=labels[ispec],kwargs...)
            clear=false
        end
    end
end

function GridVisualize.scalarplot!(vis,sys::System, sol::AbstractVector; kwargs...)
    GridVisualize.scalarplot!(vis,sys.grid,sol; kwargs...)
end
