module TestGeomspace
using Printf
using Test
using VoronoiFVM


function main()
    
    X1=VoronoiFVM.geomspace(2,3,0.2,0.2)
    X2=collect(2:0.2:3)
    @test X1 ≈  X2
    X1=VoronoiFVM.geomspace(2,3,0.1,0.2)
    X2=VoronoiFVM.geomspace(2,3,0.2,0.1)
    

    X1=VoronoiFVM.geomspace(2,3,1.0e-5,0.2)
    X2=VoronoiFVM.geomspace(2,3,0.2,1.0e-5)
    @test (X1[2]-X1[1])  ≈ (X2[end]-X2[end-1])
    @test (X2[2]-X2[1])  ≈ (X1[end]-X1[end-1])
end
end
