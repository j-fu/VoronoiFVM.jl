################################################################
"""
````
integrate(system,F,U; boundary=false)    
````

Integrate node function (same signature as reaction or storage)
 `F` of  solution vector region-wise over domain or boundary.
The result is  `nspec x nregion` vector.
"""
function integrate(system::AbstractSystem{Tv,Ti,Tm},F::Function,U::AbstractMatrix{Tu}; boundary=false) where {Tu,Tv,Ti,Tm}
    grid=system.grid
    data=system.physics.data
    nspecies=num_species(system)
    res=zeros(Tu,nspecies)
    node=Node{Tv,Ti}(system)
    nodeparams=(node,)
    bnode=BNode{Tv,Ti}(system)
    nodeparams=(bnode,)
    if isdata(data)
        nodeparams=(node,data,)
    end

    if boundary
        
        geom=grid[BFaceGeometries][1]
        bfaceregions=grid[BFaceRegions]
        nbfaceregions=maximum(bfaceregions)
        integral=zeros(Tu,nspecies,nbfaceregions)
        
        
        for ibface=1:num_bfaces(grid)
            for inode=1:num_nodes(geom)
                _fill!(bnode,inode,ibface)
                res.=zero(Tv)
                @views F(res,U[:,bnode.index],nodeparams...)
                for ispec=1:nspecies
                    if system.node_dof[ispec,bnode.index]==ispec
                        integral[ispec,bnode.region]+=system.bfacenodefactors[inode,ibface]*res[ispec]
                    end
                end
            end
        end
    else
        geom=grid[CellGeometries][1]
        cellnodes=grid[CellNodes]
        cellregions=grid[CellRegions]
        ncellregions=maximum(cellregions)
        integral=zeros(Tu,nspecies,ncellregions)
        
        
        for icell=1:num_cells(grid)
            for inode=1:num_nodes(geom)
                _fill!(node,inode,icell)
                res.=zero(Tv)
                @views F(res,U[:,node.index],nodeparams...)
                for ispec=1:nspecies
                    if system.node_dof[ispec,node.index]==ispec
                        integral[ispec,node.region]+=system.cellnodefactors[inode,icell]*res[ispec]
                    end
                end
            end
        end
    end
    
    return integral
end

