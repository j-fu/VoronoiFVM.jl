isunknownsof(u::Any, sys::AbstractSystem)=false
isunknownsof(u::DenseSolutionArray, sys::DenseSystem) = size(u) == size(sys.node_dof)
isunknownsof(u::SparseSolutionArray, sys::SparseSystem) = size(u) == size(sys.node_dof)




struct BoundaryCondition
    ispec::Int64
    ibreg::Int64
    bcfac::Float64
    bcval::Float64
end


function boundary_conditions!(sys,bcond)
    for bc ∈ bcond
        sys.boundary_factors[bc.ispec,bc.ibreg]=bc.bcfac
        sys.boundary_values[bc.ispec,bc.ibreg]=bc.bcval
    end
end


function system(grid::ExtendableGrid;
                num_species=0,
                flux::Tflux=nothing,
                reaction::Treaction=nothing,
                bc=[],
                kwargs...
              ) where {Tflux, Treaction}
    myflux= isnothing(flux) ? VoronoiFVM.nofunc : (y,u,edge)-> (y.=flux(u,edge))
    myreaction= isnothing(reaction) ? VoronoiFVM.nofunc : (y,u,node)-> (y.=reaction(u,node))
    physics=Physics(flux=myflux, reaction =myreaction)
    sys=System(grid,physics; kwargs...)
    allregions=collect(1:num_cellregions(grid))
    for ispec=1:num_species
        enable_species!(sys,ispec,allregions)
    end
    boundary_conditions!(sys,bc)
    sys
end


function VoronoiFVM.solve(sys::VoronoiFVM.System; inival=0, kwargs...)
    if isa(inival,Number)
        inival=unknowns(sys,inival=inival)
    elseif  !VoronoiFVM.isunknownsof(inival,sys)
        @error "wrong type of inival: $(typeof(inival))"
    end

    control=VoronoiFVM.NewtonControl()
    for k ∈ kwargs
        try
            setproperty!(control,k[1],k[2])
        catch e
            msg="solve: now such keyword: $(k[1])"
            @error msg
            println("Possible keywords & default values:")
            println(control)
            rethrow(e)
        end
    end
    solve(inival,sys; control=control)
end


NeumannBC(;species=nothing, boundary=nothing, value=0)=BoundaryCondition(species,boundary,0,value)
RobinBC(;species=nothing, boundary=nothing, factor=0,value=0)=BoundaryCondition(species,boundary,factor,value)
DirichletBC(;species=nothing, boundary=nothing, value=0)=BoundaryCondition(species,boundary,VoronoiFVM.Dirichlet,value)

