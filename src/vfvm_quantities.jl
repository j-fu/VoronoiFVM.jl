
abstract type AbstractQuantity{Ti<:Integer} end

struct ContinuousQuantity{Ti} <: AbstractQuantity{Ti}
    ispec::Ti
end

struct InterfaceQuantity{Ti}<: AbstractQuantity{Ti}
    ispec::Ti
end

struct DiscontinuousQuantity{Ti} <: AbstractQuantity{Ti}
    regionspec::Vector{Ti}
end


function DiscontinuousQuantity(sys::AbstractSystem{Tv,Ti,Tm},regions) where {Tv,Ti,Tm}
    nspec=num_species(sys)
    quantity=DiscontinuousQuantity{Ti}(zeros(Ti,num_cellregions(sys.grid)))
    for ireg ∈ regions
        nspec=nspec+1
        enable_species!(sys,nspec,[ireg])
        quantity.regionspec[ireg]=nspec
    end
    quantity
end

function ContinuousQuantity(sys::AbstractSystem{Tv,Ti,Tm},regions) where {Tv,Ti,Tm}
    nspec=num_species(sys)
    nspec=nspec+1
    enable_species!(sys,nspec,regions)
    ContinuousQuantity{Ti}(nspec)
end

function InterfaceQuantity(sys::AbstractSystem{Tv,Ti,Tm},regions) where {Tv,Ti,Tm}
    nspec=num_species(sys)
    nspec=nspec+1
    enable_boundary_species!(sys,nspec,regions)
    InterfaceQuantity{Ti}(nspec)
end


function subgrids(quantity::DiscontinuousQuantity, sys)
    grid=sys.grid
    subgrids=Vector{ExtendableGrid}(undef,0)
    for ireg=1:num_cellregions(grid)
        ispec=quantity.regionspec[ireg]
        if ispec>0
            push!(subgrids,subgrid(grid,[ireg]))
        end
    end
    subgrids
end

function views(U, quantity::DiscontinuousQuantity, subgrids,sys)
    grid=sys.grid
    projections=Vector[]
    j=1
    for ireg=1:num_cellregions(grid)
        ispec=quantity.regionspec[ireg]
        if ispec>0
            push!(projections,view(U[ispec,:],subgrids[j]))
            j=j+1
        end
    end
    projections
end

function views(U, quantity::ContinuousQuantity, subgrids,sys)
    grid=sys.grid
    projections=Vector[]
    j=1
    for ireg=1:num_cellregions(grid)
        push!(projections,view(U[quantity.ispec,:],subgrids[j]))
        j=j+1
    end
    projections
end



# just return the first which comes into mind.
# we need to ensure homgeneity of bregion-region structure.
# if that is true, this works.
function get_cellregion(sys,ibc)
    grid=sys.grid
    bfregions=grid[BFaceRegions]
    cregions=grid[CellRegions]
    for ibface=1:num_bfaces(sys.grid)
        if bfregions[ibface]==ibc
            bfcells=grid[BFaceCells]
            return cregions[bfcells[ibface,1]]
        end
    end
    return 0
end

function boundary_dirichlet!(sys::AbstractSystem, q::DiscontinuousQuantity, ibc, v)
    boundary_dirichlet!(sys,
                        q.regionspec[get_cellregion(sys,ibc)],
                        ibc,
                        v)
end

function boundary_dirichlet!(sys::AbstractSystem, q::ContinuousQuantity, ibc, v)
    boundary_dirichlet!(sys,
                        q.ispec,
                        ibc,
                        v)
end
