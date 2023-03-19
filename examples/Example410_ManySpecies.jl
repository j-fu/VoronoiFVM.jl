#=

# 410: Many Species
 ([source code](SOURCE_URL))

Test stationary diffusion for 50 species.

=#
module Example410_ManySpecies
using Printf
using VoronoiFVM
using SparseArrays
using ExtendableGrids
using GridVisualize
using LinearAlgebra

function main(; n = 11, nspec = 50, Plotter = nothing, unknown_storage = :dense)
    grid = simplexgrid(range(0, 1; length = n))

    function flux(f, u, edge)
        for ispec = 1:nspec
            f[ispec] = u[ispec, 1] - u[ispec, 2]
        end
    end
    physics = VoronoiFVM.Physics(; flux = flux)
    sys = VoronoiFVM.System(grid, physics; unknown_storage = unknown_storage)
    for ispec = 1:nspec
        enable_species!(sys, ispec, [1])
        boundary_dirichlet!(sys, ispec, 1, 0)
        boundary_dirichlet!(sys, ispec, 2, 1)
    end
    sol = solve(sys)
    norm(sol)
end

function test()
    testval = 13.874436925511608
    main(; unknown_storage = :sparse) ≈ testval &&
        main(; unknown_storage = :dense) ≈ testval
end

end
