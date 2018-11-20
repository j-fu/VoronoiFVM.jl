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
        "alldocs.md",
        "allindex.md",
        "Examples" => [
            "examples/nlpoisson-1spec.md",
            "examples/nlpoisson-2spec.md",
            "examples/iliq.md",
            "examples/test2d-brea-bspec.md",
            "examples/test2d-brea.md",
            "examples/test2d.md",
            "examples/test2d-rea.md"
        ]
    ]
)

