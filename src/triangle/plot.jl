#
# Plotting routines.
#

"""
$(TYPEDSIGNATURES)

Heuristic check if Plotter is PyPlot
"""
ispyplot(Plotter::Module)=isdefined(Plotter,:Gcf)

"""
$(TYPEDSIGNATURES)

Heuristic check if Plotter is Plots
"""
isplots(Plotter::Module)=isdefined(Plotter,:gr)

# Plot color scale for grid colors.
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

Plot contents of triangulateio structure.
The plot module (currently only PyPlot is possible,Plots may follow soon)
is passed as a parameter. This allows to keep the package free of heavy
plot package dependencies.
"""
function plot(Plotter::Module,tio::TriangulateIO)
    
    if ispyplot(Plotter)
        ax = Plotter.matplotlib.pyplot.gca()
        if numberofpoints(tio)==0
            return
        end
        x=tio.pointlist[1,:]
        y=tio.pointlist[2,:]
        Plotter.scatter(x,y)
        if numberoftriangles(tio)>0
            Plotter.triplot(x,y,transpose(tio.trianglelist.-1),color="k",linewidth=0.5)
        end
        if numberofsegments(tio)>0
            lines=Any[]
            rgb=Any[]
            markermax=maximum(tio.segmentmarkerlist)
            # see https://gist.github.com/gizmaa/7214002
            for i=1:numberofsegments(tio)
                x1=tio.pointlist[1,tio.segmentlist[1,i]] 
                y1=tio.pointlist[2,tio.segmentlist[1,i]]
                x2=tio.pointlist[1,tio.segmentlist[2,i]]                 
                y2=tio.pointlist[2,tio.segmentlist[2,i]]
                push!(lines,collect(zip([x1,x2],[y1,y2])))
                push!(rgb,frgb(Plotter,tio.segmentmarkerlist[i],markermax))
            end
            ax.add_collection(Plotter.matplotlib.collections.LineCollection(lines,colors=rgb,linewidth=3))
        end
    end
    
    ### Currently not working
    if isplots(Plotter)
        #=
        if p==nothing
            p=Plotter.plot()
        end
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
        if show
            Plotter.gui(p)
        end
        return p
        =#
    end
end
