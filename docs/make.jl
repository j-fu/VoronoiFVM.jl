push!(LOAD_PATH,"../src/")
using Documenter, VoronoiFVM,Markdown
makedocs(
    sitename="VoronoiFVM.jl",
    modules = [VoronoiFVM],
    clean = true,
    doctest = true,
    authors = "J. Fuhrmann",
    repo="https://github.com/j-fu/VoronoiFVM.jl",
    pages=[ 
        "Home"=>"index.md",
        "changes.md",
        "API Documentation" => [
            "grid.md",
            "physics.md",
            "system.md",
            "allindex.md",
        ],
        "Examples" => [
            "examples/Laplace.md",
            "examples/OneSpeciesNonlinearPoisson.md",
            "examples/TwoSpeciesNonlinearPoisson.md",
            "examples/IonicLiquid.md",
            "examples/NonlinearPoisson2D.md",
            "examples/NonlinearPoisson2D_Reaction.md",
            "examples/ThreeRegions1D.md",
            "examples/NonlinearPoisson2D_BoundaryReaction.md",
            "examples/NonlinearPoisson1D_BoundarySpecies.md",
            "examples/NonlinearPoisson2D_BoundarySpecies.md"
        ],
    ]
)

deploydocs(repo = "github.com/j-fu/VoronoiFVM.jl.git")


