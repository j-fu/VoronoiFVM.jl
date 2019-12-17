"""
$(TYPEDSIGNATURES)

Check if plotter is PyPlot
"""
ispyplot(Plotter)= (typeof(Plotter)==Module)&&isdefined(Plotter,:Gcf)


"""
$(TYPEDSIGNATURES)

Check if plotter is Plots
"""
isplots(Plotter)= (typeof(Plotter)==Module) && isdefined(Plotter,:gr)


"""
$(TYPEDSIGNATURES)

Check if plotter is a plotter at all.
"""
isplotter(Plotter)=isplots(Plotter)||ispyplot(Plotter)

"""
$(TYPEDSIGNATURES)

Plot color scale for grid colors.
"""
function frgb(Plotter,i,max;pastel=false)
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
    if ispyplot(Plotter)
        return (r,g,b)
    end
    if isplots(Plotter)
        return Plotter.RGB(r,g,b)
    end
end


"""
$(TYPEDSIGNATURES)

Plot grid.
"""
function plot(Plotter,grid::VoronoiFVM.Grid;
              aspect=1,
              clear=true,
              show=true,
              p=nothing)

    if ispyplot(Plotter)
        if clear
            Plotter.clf()
        end
        if dim_space(grid)==1
            xmin=minimum(grid.coord)
            xmax=maximum(grid.coord)
            h=(xmax-xmin)/20.0
            
            for icell=1:num_cells(grid)
                rgb=frgb(Plotter,grid.cellregions[icell],num_cellregions(grid))
                coord1=nodecoord(grid,cellnode(grid,1,icell))
                coord2=nodecoord(grid,cellnode(grid,2,icell))
                x1=coord1[1]
                x2=coord2[1]
                Plotter.plot([x1,x1],[-h,h],linewidth=0.5,color="k")
                Plotter.plot([x2,x2],[-h,h],linewidth=0.5,color="k")
                Plotter.plot([x1,x2],[0,0],linewidth=3.0,color=rgb)
            end
            
            for ibface=1:num_bfaces(grid)
                rgb=frgb(Plotter,grid.bfaceregions[ibface],num_bfaceregions(grid))
                coord1=nodecoord(grid,bfacenode(grid,1,ibface))
                x1=coord1[1]
                Plotter.plot([x1,x1],[-2*h,2*h],linewidth=3.0,color=rgb)
            end
        end 

        if dim_space(grid)==2
            ax=Plotter.matplotlib.pyplot.gca()
            ax.set_aspect(aspect)
            Plotter.tripcolor(tridata(grid)...,facecolors=grid.cellregions,cmap="Pastel2")
            Plotter.triplot(tridata(grid)...,color="k",linewidth=0.5)
            
            # see https://gist.github.com/gizmaa/7214002
            xc=[nodecoord(grid,bfacenode(grid,1,i)) for i=1:num_bfaces(grid)]
            yc=[nodecoord(grid,bfacenode(grid,2,i)) for i=1:num_bfaces(grid)]
            rgb=[frgb(Plotter,grid.bfaceregions[i],num_bfaceregions(grid)) for i=1:num_bfaces(grid)]
            ax.add_collection(Plotter.matplotlib.collections.LineCollection(collect(zip(xc,yc)),colors=rgb,linewidth=3))
        end
    end

    if isplots(Plotter)
        if p==nothing
            p=Plotter.plot()
        end
        if dim_space(grid)==1
            xmin=minimum(grid.coord)
            xmax=maximum(grid.coord)
            h=(xmax-xmin)/20.0
            
            for icell=1:num_cells(grid)
                rgb=frgb(Plotter,grid.cellregions[icell],num_cellregions(grid))
                coord1=nodecoord(grid,cellnode(grid,1,icell))
                coord2=nodecoord(grid,cellnode(grid,2,icell))
                x1=coord1[1]
                x2=coord2[1]
                Plotter.plot!(p,[x1,x1],[-h,h],linewidth=0.5,color=:black,label="")
                Plotter.plot!(p,[x2,x2],[-h,h],linewidth=0.5,color=:black,label="")
                Plotter.plot!(p,[x1,x2],[0,0],linewidth=3.0,color=rgb,label="")
            end
            
            for ibface=1:num_bfaces(grid)
                if grid.bfaceregions[ibface]>0
                rgb=frgb(Plotter,grid.bfaceregions[ibface],num_bfaceregions(grid))
                    coord1=nodecoord(grid,bfacenode(grid,1,ibface))
                    x1=coord1[1]
                    Plotter.plot!(p,[x1,x1],[-2*h,2*h],linewidth=3.0,color=rgb,label="")
                end
            end
        end
        
        if dim_space(grid)==2
            for icell=1:num_cells(grid)
                rgb=frgb(Plotter,grid.cellregions[icell],num_cellregions(grid),pastel=true)
                coord1=nodecoord(grid,cellnode(grid,1,icell))
                coord2=nodecoord(grid,cellnode(grid,2,icell))
                coord3=nodecoord(grid,cellnode(grid,3,icell))
                # https://github.com/JuliaPlotter/Plotter.jl/issues/605    
                tri=Plotter.Shape([coord1[1],coord2[1], coord3[1]],[coord1[2],coord2[2],coord3[2]])
                Plotter.plot!(p,tri,color=rgb,label="")
            end
            for icell=1:num_cells(grid)
                coord1=nodecoord(grid,cellnode(grid,1,icell))
                coord2=nodecoord(grid,cellnode(grid,2,icell))
                coord3=nodecoord(grid,cellnode(grid,3,icell))
                Plotter.plot!(p, [coord1[1],coord2[1]],[coord1[2],coord2[2]]  ,linewidth=0.5,color=:black,label="")
                Plotter.plot!(p, [coord1[1],coord3[1]],[coord1[2],coord3[2]]  ,linewidth=0.5,color=:black,label="")
                Plotter.plot!(p, [coord2[1],coord3[1]],[coord2[2],coord3[2]]  ,linewidth=0.5,color=:black,label="")
            end
            for ibface=1:num_bfaces(grid)
                rgb=frgb(Plotter,grid.bfaceregions[ibface],num_bfaceregions(grid))
                coord1=nodecoord(grid,bfacenode(grid,1,ibface))
                coord2=nodecoord(grid,bfacenode(grid,2,ibface))
                Plotter.plot!(p,[coord1[1],coord2[1]],[coord1[2],coord2[2]]  ,linewidth=5,color=rgb,label="")
            end
        end
        if show
            Plotter.gui(p)
        end
        return p
    end
end


"""
$(TYPEDSIGNATURES)

Plot subgrid.
"""
function plot(Plotter,subgrid::VoronoiFVM.SubGrid;
              clear=true,
              show=true)
    
    if ispyplot(Plotter)
        if clear
            Plotter.clf()
        end
        
        if dim_space(subgrid.parent)==1
        xmin=minimum(subgrid.coord)
            xmax=maximum(subgrid.coord)
            h=(xmax-xmin)/20.0
            Plotter.xlim(xmin-h,xmax+h)
            Plotter.ylim(-20*h,20*h)
            
            for icell=1:num_cells(subgrid)
                coord1=nodecoord(subgrid,subgrid.cellnodes[1,icell])
                coord2=nodecoord(subgrid,subgrid.cellnodes[2,icell])
                x1=coord1[1]
                x2=coord2[1]
                Plotter.plot([x1,x1],[-h,h],linewidth=0.5,color="k")
                Plotter.plot([x2,x2],[-h,h],linewidth=0.5,color="k")
                Plotter.plot([x1,x2],[0,0],linewidth=3.0,color="k")
            end
        end
    end


    if isplots(Plotter)
        p=Plotter.plot()
        if dim_space(subgrid.parent)==1
            xmin=minimum(subgrid.coord)
            xmax=maximum(subgrid.coord)
            h=(xmax-xmin)/20.0
            Plotter.xlim!(p,xmin-h,xmax+h)
            Plotter.ylim!(p,-20*h,20*h)
            
            for icell=1:num_cells(subgrid)
                coord1=nodecoord(subgrid,subgrid.cellnodes[1,icell])
                coord2=nodecoord(subgrid,subgrid.cellnodes[2,icell])
                x1=coord1[1]
                x2=coord2[1]
                Plotter.plot!(p,[x1,x1],[-h,h],linewidth=0.5,color=:black)
                Plotter.plot!(p,[x2,x2],[-h,h],linewidth=0.5,color=:black)
                Plotter.plot!(p,[x1,x2],[0,0],linewidth=3.0,color=:black)
            end
        end
        if show
            Plotter.gui(p)
        end
        return p
    end
end

"""
$(TYPEDSIGNATURES)

Plot array as piecewise linear function on grid.
"""
function plot(Plotter,grid::VoronoiFVM.AbstractGrid, U::AbstractArray;
              color=(0,0,0),
              cmap="hot",
              label="",
              levels=10,
              aspect=1,
              clear=true,
              show=true)

    if ispyplot(Plotter)
        if clear
            Plotter.clf()
        end
        if dim_space(grid)==1
            for icell=1:num_cells(grid)
                i1=grid.cellnodes[1,icell]
                i2=grid.cellnodes[2,icell]
                x1=grid.coord[1,i1]
                x2=grid.coord[1,i2]
                if icell==1 && label !=""
                    Plotter.plot([x1,x2],[U[i1],U[i2]],color=color,label=label)
                else
                    Plotter.plot([x1,x2],[U[i1],U[i2]],color=color)
                end                
            end
        end
        
        if dim_space(grid)==2
            ax=Plotter.matplotlib.pyplot.gca()
            ax.set_aspect(aspect)
            plotted=Plotter.tricontourf(tridata(grid)...,U;levels=levels,cmap=cmap)
            Plotter.tricontour(tridata(grid)...,U,colors="k",levels=levels)
            return plotted
        end
    end

    if isplots(Plotter)
        p=Plotter.plot()
        if dim_space(grid)==1
            for icell=1:num_cells(grid)
                i1=grid.cellnodes[1,icell]
                i2=grid.cellnodes[2,icell]
                x1=grid.coord[1,i1]
                x2=grid.coord[1,i2]
                if icell==1
                    Plotter.plot!(p,[x1,x2],[U[i1],U[i2]],color=Plotter.RGB(color...),label=label)
                else
                    Plotter.plot!(p,[x1,x2],[U[i1],U[i2]],color=Plotter.RGB(color...),label="")
                end                
            end
        end
        if dim_space(grid)==2
            println("Not available for Plots, see e.g. https://github.com/JuliaPlots/Plots.jl/issues/392")
        end
        if show
            Plotter.gui(p)
        end
        return p
    end
end

    

"""
$(TYPEDSIGNATURES)

Plot array as piecewise constant function on subgrid
"""
function plot(Plotter,grid::VoronoiFVM.SubGrid, U::AbstractArray;
              color=(0,0,0),
              label="",
              clear=true,
              show=true,
              p=nothing)

    if ispyplot(Plotter)
        if clear
            Plotter.clf()
        end
        if dim_space(grid)==1
            for icell=1:num_cells(grid)
                i1=grid.cellnodes[1,icell]
                i2=grid.cellnodes[2,icell]
                x1=grid.coord[1,i1]
                x2=grid.coord[1,i2]
                ip1=grid.node_in_parent[i1]
                ip2=grid.node_in_parent[i2]
                
                if icell==1 && label !=""
                    Plotter.plot([x1,x2],[U[ip1],U[ip2]],color=color,label=label)
                else
                    Plotter.plot([x1,x2],[U[ip1],U[ip2]],color=color)
                end                
            end
        end
    end

    if isplots(Plotter)
        if p==nothing
            p=Plot()
        end
        if dim_space(grid)==1
            for icell=1:num_cells(grid)
                i1=grid.cellnodes[1,icell]
                i2=grid.cellnodes[2,icell]
                x1=grid.coord[1,i1]
                x2=grid.coord[1,i2]
                ip1=grid.node_in_parent[i1]
                ip2=grid.node_in_parent[i2]
                
                if icell==1 && label !=""
                    Plotter.plot!(p,[x1,x2],[U[ip1],U[ip2]],label=label,linecolor=Plotter.RGB(color...))
                else
                    Plotter.plot!(p,[x1,x2],[U[ip1],U[ip2]],label="",linecolor=Plotter.RGB(color...))
                end                
            end
        end
        if show
            Plotter.gui(p)
        end
        return p
    end
end

