import PyPlot

function frgb(i,max;pastel=false)
    x=Float64(i-1)/Float64(max)
    if (x<0.5)
        r=1.0-2.0*x
        g=2.0*x
        b=0.0
    else
        r=0.0
        g=2.0-2.0*x
        b=2.0*x-1.0
    end
    if pastel
        r=0.5+0.5*r
        g=0.5+0.5*g
        b=0.5+0.5*b
    end
    return (r,g,b)
end



function fvmplot(grid::Grid)
    
    
    if dim_space(grid)==1
        xmin=minimum(grid.coord)
        xmax=maximum(grid.coord)
        h=(xmax-xmin)/20.0

        for icell=1:num_cells(grid)
            rgb=frgb(grid.cellregions[icell],num_cellregions(grid))
            coord1=nodecoord(grid,cellnode(grid,1,icell))
            coord2=nodecoord(grid,cellnode(grid,2,icell))
            x1=coord1[1]
            x2=coord2[1]
            PyPlot.plot([x1,x1],[-h,h],linewidth=0.5,color="k")
            PyPlot.plot([x2,x2],[-h,h],linewidth=0.5,color="k")
            PyPlot.plot([x1,x2],[0,0],linewidth=3.0,color=rgb)
        end
        
        for ibface=1:num_bfaces(grid)
            rgb=frgb(grid.bfaceregions[ibface],num_bfaceregions(grid))
            coord1=nodecoord(grid,bfacenode(grid,1,ibface))
            x1=coord1[1]
            PyPlot.plot([x1,x1],[-2*h,2*h],linewidth=3.0,color=rgb)
        end
    end


    if dim_space(grid)==2
        for icell=1:num_cells(grid)
            rgb=frgb(grid.cellregions[icell],num_cellregions(grid),pastel=true)
            coord1=nodecoord(grid,cellnode(grid,1,icell))
            coord2=nodecoord(grid,cellnode(grid,2,icell))
            coord3=nodecoord(grid,cellnode(grid,3,icell))
            PyPlot.fill( [coord1[1],coord2[1], coord3[1]],[coord1[2],coord2[2],coord3[2]]  ,linewidth=0.5,color=rgb)
        end
        for icell=1:num_cells(grid)
            coord1=nodecoord(grid,cellnode(grid,1,icell))
            coord2=nodecoord(grid,cellnode(grid,2,icell))
            coord3=nodecoord(grid,cellnode(grid,3,icell))
            PyPlot.plot( [coord1[1],coord2[1]],[coord1[2],coord2[2]]  ,linewidth=0.5,color="k")
            PyPlot.plot( [coord1[1],coord3[1]],[coord1[2],coord3[2]]  ,linewidth=0.5,color="k")
            PyPlot.plot( [coord2[1],coord3[1]],[coord2[2],coord3[2]]  ,linewidth=0.5,color="k")
        end
        for ibface=1:num_bfaces(grid)
            rgb=frgb(grid.bfaceregions[ibface],num_bfaceregions(grid))
            coord1=nodecoord(grid,bfacenode(grid,1,ibface))
            coord2=nodecoord(grid,bfacenode(grid,2,ibface))
            PyPlot.plot( [coord1[1],coord2[1]],[coord1[2],coord2[2]]  ,linewidth=5,color=rgb)
        end

    end

end


function fvmplot(subgrid::SubGrid)
    if dim_space(subgrid.parent)==1
        xmin=minimum(subgrid.coord)
        xmax=maximum(subgrid.coord)
        h=(xmax-xmin)/20.0
        PyPlot.xlim(xmin-h,xmax+h)
        PyPlot.ylim(-20*h,20*h)
        
        for icell=1:num_cells(subgrid)
            coord1=nodecoord(subgrid,subgrid.cellnodes[1,icell])
            coord2=nodecoord(subgrid,subgrid.cellnodes[2,icell])
            x1=coord1[1]
            x2=coord2[1]
            PyPlot.plot([x1,x1],[-h,h],linewidth=0.5,color="k")
            PyPlot.plot([x2,x2],[-h,h],linewidth=0.5,color="k")
            PyPlot.plot([x1,x2],[0,0],linewidth=3.0,color="k")
        end
    end
end



function fvmplot(grid::AbstractGrid, U::AbstractArray; color=(0,0,0),label="")
    if dim_space(grid)==1
        for icell=1:num_cells(grid)
            i1=grid.cellnodes[1,icell]
            i2=grid.cellnodes[2,icell]
            x1=grid.coord[1,i1]
            x2=grid.coord[1,i2]
            if icell==1 && label !=""
                PyPlot.plot([x1,x2],[U[i1],U[i2]],color=color,label=label)
            else
                PyPlot.plot([x1,x2],[U[i1],U[i2]],color=color)
            end                
        end
    end
end


function fvmplot(grid::SubGrid, U::Array{Tv,1}; color=(0,0,0),label="") where Tv
    if dim_space(grid)==1
        for icell=1:num_cells(grid)
            i1=grid.cellnodes[1,icell]
            i2=grid.cellnodes[2,icell]
            x1=grid.coord[1,i1]
            x2=grid.coord[1,i2]
            ip1=grid.node_in_parent[i1]
            ip2=grid.node_in_parent[i2]

            if icell==1 && label !=""
                PyPlot.plot([x1,x2],[U[ip1],U[ip2]],color=color,label=label)
            else
                PyPlot.plot([x1,x2],[U[ip1],U[ip2]],color=color)
            end                
        end
    end
end


