import PyPlot
import Colors


function fvmplot(grid::Grid)
    rgb_sequence(c) = (Colors.red(c), Colors.green(c), Colors.blue(c))
    ccols=Colors.distinguishable_colors(ncellregions(grid),lchoices=50:50,cchoices=50:50)
    bcols=Colors.distinguishable_colors(nbfaceregions(grid),lchoices=50:50,cchoices=50:50)
    
    if spacedim(grid)==1
        xmin=minimum(grid.nodecoord)
        xmax=maximum(grid.nodecoord)
        h=(xmax-xmin)/20.0

        for icell=1:ncells(grid)
            rgb=rgb_sequence(ccols[grid.cellregions[icell]])
            coord1=nodecoord(grid,cellnodes(grid,1,icell))
            coord2=nodecoord(grid,cellnodes(grid,2,icell))
            x1=coord1[1]
            x2=coord2[1]
            PyPlot.plot([x1,x1],[-h,h],linewidth=0.5,color="k")
            PyPlot.plot([x2,x2],[-h,h],linewidth=0.5,color="k")
            PyPlot.plot([x1,x2],[0,0],linewidth=3.0,color=rgb)
        end
        
        for ibface=1:nbfaces(grid)
            rgb=rgb_sequence(bcols[grid.bfaceregions[ibface]])
            coord1=nodecoord(grid,bfacenodes(grid,1,ibface))
            x1=coord1[1]
            PyPlot.plot([x1,x1],[-2*h,2*h],linewidth=3.0,color=rgb)
        end
    end
end


function fvmplot(subgrid::SubGrid)
    if spacedim(subgrid.parent)==1
        xmin=minimum(subgrid.nodecoord)
        xmax=maximum(subgrid.nodecoord)
        h=(xmax-xmin)/20.0
        PyPlot.xlim(xmin-h,xmax+h)
        PyPlot.ylim(-20*h,20*h)
        
        for icell=1:ncells(subgrid)
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
    if spacedim(grid)==1
        for icell=1:ncells(grid)
            i1=grid.cellnodes[1,icell]
            i2=grid.cellnodes[2,icell]
            x1=grid.nodecoord[1,i1]
            x2=grid.nodecoord[1,i2]
            if icell==1 && label !=""
                PyPlot.plot([x1,x2],[U[i1],U[i2]],color=color,label=label)
            else
                PyPlot.plot([x1,x2],[U[i1],U[i2]],color=color)
            end                
        end
    end
end
