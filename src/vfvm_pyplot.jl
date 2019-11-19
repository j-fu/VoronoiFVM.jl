# This is "dead code" meant for some ideas.
import PyPlot



function pyplot_frgb(i,max;pastel=false)
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



function fvmpyplot(grid::Grid; aspect=1)
    if dim_space(grid)==1
        xmin=minimum(grid.coord)
        xmax=maximum(grid.coord)
        h=(xmax-xmin)/20.0

        for icell=1:num_cells(grid)
            rgb=pyplot_frgb(grid.cellregions[icell],num_cellregions(grid))
            coord1=nodecoord(grid,cellnode(grid,1,icell))
            coord2=nodecoord(grid,cellnode(grid,2,icell))
            x1=coord1[1]
            x2=coord2[1]
            PyPlot.plot([x1,x1],[-h,h],linewidth=0.5,color="k")
            PyPlot.plot([x2,x2],[-h,h],linewidth=0.5,color="k")
            PyPlot.plot([x1,x2],[0,0],linewidth=3.0,color=rgb)
        end
        
        for ibface=1:num_bfaces(grid)
            rgb=pylot_frgb(grid.bfaceregions[ibface],num_bfaceregions(grid))
            coord1=nodecoord(grid,bfacenode(grid,1,ibface))
            x1=coord1[1]
            PyPlot.plot([x1,x1],[-2*h,2*h],linewidth=3.0,color=rgb)
        end
    end 

    if dim_space(grid)==2
        PyPlot.axes(aspect=aspect)
        PyPlot.tripcolor(tridata(grid)...,facecolors=grid.cellregions,cmap="Pastel2")
        PyPlot.triplot(tridata(grid)...,color="k",linewidth=0.5)
        for ibface=1:num_bfaces(grid)
            rgb=pyplot_frgb(grid.bfaceregions[ibface],num_bfaceregions(grid))
            coord1=nodecoord(grid,bfacenode(grid,1,ibface))
            coord2=nodecoord(grid,bfacenode(grid,2,ibface))
            PyPlot.plot( [coord1[1],coord2[1]],[coord1[2],coord2[2]]  ,linewidth=3,color=rgb)
        end
    end
end


function fvmpyplot(subgrid::SubGrid)
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



function fvmpyplot(grid::AbstractGrid, U::AbstractArray; color=(0,0,0),cmap="hot",label="",levels=10,aspect=1)
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
    if dim_space(grid)==2
        PyPlot.axes(aspect=aspect)
        PyPlot.tricontourf(tridata(grid)...,U;levels=levels,cmap=cmap)
        PyPlot.tricontour(tridata(grid)...,U,colors="k",levels=levels)
    end
end

    


function fvmpyplot(grid::SubGrid, U::Array{Tv,1}; color=(0,0,0),label="") where Tv
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


