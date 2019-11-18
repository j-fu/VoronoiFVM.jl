module test_gridconstruct
using VoronoiFVM

function test()
    grid=VoronoiFVM.Grid(points=[1 0 ; 0 1 ; -1 0 ; 0 -1],
                         bfaces=[1 2 ; 2 3 ; 3 4 ; 4 1 ],
                         bfaceregions=[1, 2, 3, 4],
                         regionpoints=[0.5 0.5;],
                         regionnumbers=[1],
                         regionvolumes=[0.01],
                         flags="pAaqQD")
    num_nodes(grid)==183 && num_cells(grid)==308 && num_bfaces(grid)==56
end
end
