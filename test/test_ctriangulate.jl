module test_ctriangulate
using VoronoiFVM


function test()
    nodes=Matrix{Cdouble}([1.0 0.0 ; 0.0 1.0 ; -1.0 0.0 ; 0.0 -1.0]')
    faces=Matrix{Cint}([1 2 ; 2 3 ; 3 4 ; 4 1 ]')
    faceregions=Matrix{Cint}([1 2 3 4]')
    regionpoints=Matrix{Cdouble}([0.5 0.5 1 0.01;]')
    regionnumbers=[1]
    triin=VoronoiFVM.Triangle.CTriangulateIO()
    triout=VoronoiFVM.Triangle.CTriangulateIO()
    vorout=VoronoiFVM.Triangle.CTriangulateIO()
    triin.numberofpoints=Cint(size(nodes,2))
    triin.pointlist=pointer(nodes)
    triin.numberofsegments=size(faces,2)
    triin.segmentlist=pointer(faces)
    triin.segmentmarkerlist=pointer(faceregions)
    triin.numberofregions=size(regionpoints,2)
    triin.regionlist=pointer(regionpoints)
    
    VoronoiFVM.Triangle.triangulate("paAqQ",triin,triout,vorout)
    points = convert(Array{Float64,2}, Base.unsafe_wrap(Array, triout.pointlist, (2,Int(triout.numberofpoints)), own=true))
    cells  = convert(Array{Int32,2}, Base.unsafe_wrap(Array, triout.trianglelist, (2,Int(triout.numberoftriangles)), own=true))
    bfaces = convert(Array{Int32,2}, Base.unsafe_wrap(Array, triout.segmentlist, (2,Int(triout.numberofsegments)), own=true))
    cellregions=convert(Array{Float64,1}, Base.unsafe_wrap(Array, triout.triangleattributelist, (Int(triout.numberoftriangles)), own=true))
    bfaceregions=convert(Array{Int32,1}, Base.unsafe_wrap(Array, triout.segmentmarkerlist, (Int(triout.numberofsegments)), own=true))
    cellregions=Vector{Int32}(cellregions)
    
    grid=VoronoiFVM.Grid(points,cells,cellregions,bfaces,bfaceregions)
    
    num_nodes(grid)==177 && num_cells(grid)==319 && num_bfaces(grid)==33
end
end
