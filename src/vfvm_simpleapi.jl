isunknownsof(u::Any, sys::AbstractSystem)=false
isunknownsof(u::DenseSolutionArray, sys::DenseSystem) = size(u) == size(sys.node_dof)
isunknownsof(u::SparseSolutionArray, sys::SparseSystem) = size(u) == size(sys.node_dof)



boundary_dirichlet!(y,u,bnode,ispec,ireg,val) =      bnode.region == ireg ? y[ispec] += Dirichlet*(u[ispec]-val) : nothing
boundary_robin!(y,u,bnode,ispec,ireg,fac,val) =      bnode.region == ireg ? y[ispec] += fac*u[ispec]-val : nothing
boundary_neumann!(y,u,bnode,ispec,ireg,fac,val) =    bnode.region == ireg ? y[ispec] -= val : nothing

boundary_dirichlet!(y,u,bnode; species=1, region=1,value=0) =      bnode.region == region ? y[species] += Dirichlet*(u[species]-value) : nothing


function system(grid::ExtendableGrid;
                kwargs...
                )
    physics=Physics(; kwargs...)
    sys=System(grid,physics; kwargs...)
    allregions=collect(1:num_cellregions(grid))
    num_species=kwargs[:num_species]
    for ispec=1:num_species
        enable_species!(sys,ispec,allregions)
    end
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

