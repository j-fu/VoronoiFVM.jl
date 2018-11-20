push!(LOAD_PATH,"../src/")
using Documenter, TwoPointFluxFVM
makedocs(
    sitename="TwoPointFluxFVM.jl",
    modules = [TwoPointFluxFVM],
    clean = true,
    authors = "J. Fuhrmann",
    pages=["Home"=>"index.md",
           "alldocs.md",
           "allindex.md"
           ]
)

