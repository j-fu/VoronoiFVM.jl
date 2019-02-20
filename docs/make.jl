push!(LOAD_PATH,"../src/")
using Documenter, TwoPointFluxFVM
makedocs(
    sitename="TwoPointFluxFVM.jl",
    modules = [TwoPointFluxFVM],
    clean = true,
    authors = "J. Fuhrmann",
    pages=[
        "Home"=>"index.md",
        "install.md",
        "changes.md",
        "alldocs.md",
        "allindex.md",
        "Examples" => [
            "examples/OneSpeciesNonlinearPoisson.md",
            "examples/TwoSpeciesNonlinearPoisson.md",
            "examples/IonicLiquid.md",
            "examples/NonlinearPoisson2D.md",
            "examples/NonlinearPoisson2D_Reaction.md",
            "examples/NonlinearPoisson2D_BoundaryReaction.md",
            "examples/NonlinearPoisson1D_BoundarySpecies.md",
            "examples/NonlinearPoisson2D_BoundarySpecies.md"
        ]
    ]
)

