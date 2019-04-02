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
        PyPlot.xlim(xmin-h,xmax+h)
        PyPlot.ylim(-20*h,20*h)

        for icell=1:ncells(grid)
            rgb=rgb_sequence(ccols[grid.cellregions[icell]])
            coord1=nodecoord(grid,cellnodes(grid,1,icell))
            coord2=nodecoord(grid,cellnodes(grid,2,icell))
            x1=coord1[1]
            x2=coord2[1]
            PyPlot.plot([x1,x1],[-h,h],linewidth=0.5,color="b")
            PyPlot.plot([x2,x2],[-h,h],linewidth=0.5,color="b")
            PyPlot.plot([x1,x2],[0,0],linewidth=3.0,color=rgb)
        end
        
        for ibface=1:nbfaces(grid)
            rgb=rgb_sequence(bcols[grid.bfaceregions[ibface]])
            coord1=nodecoord(grid,bfacenodes(grid,1,ibface))
            x1=coord1[1]
            println(x1)
            PyPlot.plot([x1,x1],[-2*h,2*h],linewidth=3.0,color=rgb)
        end
    end
    
end
