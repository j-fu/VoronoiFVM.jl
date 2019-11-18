module test_triangulate
using VoronoiFVM

function test()
    triin=VoronoiFVM.TriangulateIO()
    triin.pointlist=Matrix{Float64}([1.0 0.0 ; 0.0 1.0 ; -1.0 0.0 ; 0.0 -1.0]')
    triin.segmentlist=Matrix{Int32}([1 2 ; 2 3 ; 3 4 ; 4 1 ]')
    triin.segmentmarkerlist=Vector{Int32}([1, 2, 3, 4])
    triin.regionlist=Matrix{Float64}([0.5 0.5 1 0.01;]')
    grid=VoronoiFVM.Grid("paAqQ",triin)
    num_nodes(grid)==177 && num_cells(grid)==319 && num_bfaces(grid)==33
end
end
