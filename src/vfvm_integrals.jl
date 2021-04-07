################################################################
"""
````
integrate(system,F,U; boundary=false)    
````

Integrate node function (same signature as reaction or storage)
 `F` of  solution vector over domain or boundary 
The result contains the integral for each species separately.
"""
function integrate(system::AbstractSystem{Tv,Ti,Tm},F::Function,U::AbstractMatrix{Tu}; boundary=false) where {Tu,Tv,Ti,Tm}
    grid=system.grid
    data=system.physics.data
    nspecies=num_species(system)
    res=zeros(Tu,nspecies)
    node=Node{Tv,Ti}(system)
    nodeparams=(node,)
    if isdata(data)
        nodeparams=(node,data,)
    end

    

    csys=grid[CoordinateSystem]
    coord=grid[Coordinates]


    if boundary
        
        geom=grid[BFaceGeometries][1]
        bfacenodes=grid[BFaceNodes]
        bfaceregions=grid[BFaceRegions]
        nbfaceregions=maximum(bfaceregions)
        integral=zeros(Tu,nspecies,nbfaceregions)
        
        
        for ibface=1:num_bfaces(grid)
            for inode=1:num_nodes(geom)
                _fill!(node,bfacenodes,bfaceregions,inode,ibface)
                res.=zero(Tv)
                @views F(res,U[:,node.index],nodeparams...)
                for ispec=1:nspecies
                    if system.node_dof[ispec,node.index]==ispec
                        integral[ispec,node.region]+=system.bfacenodefactors[inode,ibface]*res[ispec]
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
                _fill!(node,cellnodes,cellregions,inode,icell)
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

