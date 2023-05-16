module test_norms
using VoronoiFVM
using ExtendableGrids
using Test
using LinearAlgebra

function grid(X, dim)
    if dim == 1
        simplexgrid(X)
    elseif dim == 2
        simplexgrid(X, X)
    else
        simplexgrid(X, X, X)
    end
end


function test_const(; dim = 2, c = 1.0, nrm = l2norm, assembly = :edgewise, h = 0.1)
    X = 0:h:1
    g = grid(X, dim)
    sys = VoronoiFVM.System(g; species = [1], assembly)
    VoronoiFVM._complete!(sys)
    F = map(c, sys)
    nrm(sys, F)
end


function test_lin(; dim = 2, c = 1.0, nrm = h1seminorm, assembly = :edgewise, h = 0.1)
    X = 0:h:1
    g = grid(X, dim)
    sys = VoronoiFVM.System(g; species = [1], assembly)
    VoronoiFVM._complete!(sys)
    F = unknowns(sys)
    F[1, :] .= map((x...) -> (c * sum(x) / dim), g)
    nrm(sys, F)
end

function test_transient(;
    dim = 2,
    c = 1.0,
    nrm = l2h1seminorm,
    assembly = :edgewise,
    h = 0.1,
)
    X = 0:h:1
    g = grid(X, dim)
    sys = VoronoiFVM.System(g; species = [1], assembly)
    VoronoiFVM._complete!(sys)
    F = unknowns(sys)
    F[1, :] .= map((x...) -> (c * sum(x) / dim), g)
    U = TransientSolution(X[1], F)
    for i = 2:length(X)
        append!(U, X[i], F)
    end

    nrm(sys, U)
end

function test()
    for assembly in (:edgewise, :cellwise)
        for dim = 1:3
            for c in [0.5, 1, 2.0]
                @test test_const(; nrm = l2norm, dim, c, assembly) ≈ c
                @test test_const(; nrm = h1seminorm, dim, c, assembly) ≈ 0
                @test test_lin(; nrm = h1seminorm, dim, c, assembly) ≈ c
                @test test_transient(; nrm = l2h1seminorm, dim, c, assembly) ≈ c
            end
        end
    end
    true
end

end
