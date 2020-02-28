#=
Definition of methods to edit grid region numbers
=#

######################################################
"""
$(TYPEDSIGNATURES)

Edit region numbers of grid cells via rectangular mask.
"""
function cellmask!(grid::Grid,
                   maskmin::AbstractArray,
                   maskmax::AbstractArray,
                   ireg::Int;
                   tol=1.0e-10)
    xmaskmin=maskmin.-tol
    xmaskmax=maskmax.+tol
    for icell=1:num_cells(grid)
        in_region=true
        for inode=1:num_nodes_per_cell(grid)
            coord=nodecoord(grid,cellnode(grid,inode,icell))
            for idim=1:dim_space(grid)
                if coord[idim]<xmaskmin[idim]
                    in_region=false
                elseif coord[idim]>xmaskmax[idim]
                    in_region=false
                end
            end
        end
        if in_region
            grid.cellregions[icell]=ireg
        end
    end
    grid.num_cellregions=max(num_cellregions(grid),ireg)
end


######################################################
"""
$(TYPEDSIGNATURES)

Edit region numbers of grid  boundary facets  via rectangular mask.
Currently, only for 1D grids, inner boundaries can be added.
"""
function bfacemask!(grid::Grid,
                    maskmin::AbstractArray,
                    maskmax::AbstractArray,
                    ireg::Int;
                    tol=1.0e-10)

    
    
    xmaskmin=maskmin.-tol
    xmaskmax=maskmax.+tol
    
    function isbface(ix)
        for ibface=1:num_bfaces(grid)
            if grid.bfacenodes[1,ibface]==ix
                return ibface
            end
            return 0
        end
    end
    if dim_space(grid)==1
        Ti=eltype(grid.bfacenodes)
        bfacenodes=ElasticArray{Ti,2}(grid.bfacenodes)
        for inode=1:num_nodes(grid)
            x=grid.coord[1,inode]
            if x>xmaskmin[1] && x<xmaskmax[1]
                ibface=isbface(inode)
                if ibface>0
                    grid.bfaceregions[ibface]=ireg
                else
                    ibface=length(grid.bfaceregions)+1
                    push!(grid.bfaceregions,ireg)
                    append!(bfacenodes,[inode])
                end
            end
        end
        grid.bfacenodes=Array{Ti,2}(bfacenodes)
    else
        for ibface=1:num_bfaces(grid)
            in_region=true
            for inode=1:num_nodes_per_bface(grid)
                coord=nodecoord(grid,bfacenode(grid,inode,ibface))
                for idim=1:dim_space(grid)
                    if coord[idim]<xmaskmin[idim]
                        in_region=false
                    elseif coord[idim]>xmaskmax[idim]
                        in_region=false
                    end
                end
            end
            if in_region
                grid.bfaceregions[ibface]=ireg
            end
        end
    end
        
    grid.num_bfaceregions=max(num_bfaceregions(grid),ireg)
    return grid
end
