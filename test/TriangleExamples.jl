# TriangleExamples
# ================
#
# These are a couple of more or less useful first examples which at the
# same time serve as tests. They can be loaded into Julia by
# `includet("examples/TriangleExamples.jl")` and run via e.g.
# `TriangleExamples.main(plotter=PyPlot,example="convexhull")`
# 
# 
module TriangleExamples

# Include TriangleRaw and Test module

using VoronoiFVM.Triangle
using Test


#
# Plot a pair of input and output triangulateio structs
#
function plotpair(Plotter::Module, triin, triout;voronoi=nothing,title="")
    if ispyplot(Plotter)
        PyPlot=Plotter
        PyPlot.clf()
        PyPlot.suptitle(title)
        PyPlot.subplot(121)
        PyPlot.title("In")
        Triangle.plot(PyPlot,triin)
        PyPlot.subplot(122)
        PyPlot.title("Out")
        Triangle.plot(PyPlot,triout,voronoi=voronoi)
    end
end


#
# This function can be called interactively
# with `Plotter =PyPlot` and example set to one
# of the example strings. "Plots" is in the making...
#
# At the same time it is used as part of the
# runtime tests.
function main(;Plotter=Triangle, example="all")
   
    do_example(ex)= example==ex || example=="all"

    # Handle plot setup
    if ispyplot(Plotter)
        fig = Plotter.matplotlib.pyplot.gcf()
        fig.set_size_inches(10,5)
    end
    
    # Delaunay triangulation of convex hull
    if do_example("convexhull")
        triin=Triangle.TriangulateIO()
        triin.pointlist=rand(Cdouble,2,20)
        (triout, vorout)=triangulate("c", triin)
        @test numberofpoints(triin)==numberofpoints(triout)
        @test numberofsegments(triout)>0
        plotpair(Plotter,triin,triout,title=example)
    end

    # Delaunay triangulation of convex hull
    if do_example("vconvexhull")
        triin=Triangle.TriangulateIO()
        triin.pointlist=[ 1 2 ; 3 4]'
        (triout, vorout)=triangulate("cv", triin)
        plotpair(Plotter,triin,triout,voronoi=vorout,title=example)
    end


    
    # Constrained Delaunay triangulation 
    if do_example("cdt")
        triin=Triangle.TriangulateIO()
        triin.pointlist=rand(Cdouble,2,20)
        triin.segmentlist=Matrix{Cint}([1 20; 9 10]')
        triin.segmentmarkerlist=Vector{Cint}([2,3])
        (triout, vorout)=triangulate("pc", triin)
        @test numberofpoints(triin)<=numberofpoints(triout)
        @test numberofsegments(triout)>0
        plotpair(Plotter,triin,triout,title=example)
    end
    
    # Delaunay triangulation of pointset
    if do_example("dcdt")
        triin=Triangle.TriangulateIO()
        triin.pointlist=rand(Cdouble,2,20)
        triin.segmentlist=Matrix{Cint}([1 20; 9 10]')
        triin.segmentmarkerlist=Vector{Cint}([2,3])
        (triout, vorout)=triangulate("Dpc", triin)
        @test numberofpoints(triin)<=numberofpoints(triout)
        @test numberofsegments(triout)>0
        plotpair(Plotter,triin,triout,title=example)
    end
    
    # Delaunay triangulation of domain
    if do_example("domain")
        triin=Triangle.TriangulateIO()
        triin.pointlist=Matrix{Cdouble}([0.0 0.0 ; 1.0 0.0 ; 0.9  0.9 ; 0.0 1.0]')
        triin.segmentlist=Matrix{Cint}([1 2 ; 2 3 ; 3 4 ; 4 1 ]')
        triin.segmentmarkerlist=Vector{Int32}([1, 2, 3, 4])
        (triout, vorout)=triangulate("pa0.01", triin)
        @test numberofpoints(triout)==87
        @test numberofsegments(triout)==21
        @test numberoftriangles(triout)==151
        plotpair(Plotter,triin,triout,title=example)
    end
    # Boundary conforming Delaunay triangulation of domain
    if do_example("ddomain")
        triin=Triangle.TriangulateIO()
        triin.pointlist=Matrix{Cdouble}([0.0 0.0 ; 1.0 0.0 ; 0.9  0.9 ; 0.0 1.0]')
        triin.segmentlist=Matrix{Cint}([1 2 ; 2 3 ; 3 4 ; 4 1 ]')
        triin.segmentmarkerlist=Vector{Int32}([1, 2, 3, 4])
        (triout, vorout)=triangulate("pa0.01D", triin)
        @test numberofpoints(triout)==84
        @test numberofsegments(triout)==32
        @test numberoftriangles(triout)==134
        plotpair(Plotter,triin,triout,title=example)
    end

    # Example with local refinement calling user triunsuitable function
    if do_example("localref")

        center_x=5.0
        center_y=5.0
        localdist=1.0
        localarea=0.01
        minarea=0.5

        function unsuitable(x1,y1,x2,y2,x3,y3,area)
            bary_x=(x1+x2+x3)/3.0
            bary_y=(y1+y2+y3)/3.0
            dx=bary_x-center_x
            dy=bary_y-center_y
            dist=dx^2+dy^2
            
            (dist^2<localdist^2 && area>localarea) || (area>minarea)
        end
        
        triunsuitable(unsuitable)
        triin=Triangle.TriangulateIO()
        triin.pointlist=Matrix{Cdouble}([0.0 0.0 ; 10.0 0.0 ; 10.0  10.0 ; 0.0 10.0]')
        triin.segmentlist=Matrix{Cint}([1 2 ; 2 3 ; 3 4 ; 4 1 ]')
        triin.segmentmarkerlist=Vector{Int32}([1, 2, 3, 4])
        (triout, vorout)=triangulate("puaD", triin)
        @test numberofpoints(triout)==420
        @test numberofsegments(triout)==48
        @test numberoftriangles(triout)==790
        plotpair(Plotter,triin,triout,title=example)
    end
    true
end

# 
# Called by runtest.
#
function test()
    main(;Plotter=Triangle, example="all")
end

# End of module
end
