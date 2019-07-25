using Documenter, VoronoiFVM,Literate



output_dir = joinpath(@__DIR__, "src","examples")
example_sources=readdir(joinpath(@__DIR__,"..","examples"))
for example_source in example_sources
    println("$(example_source):")
    Literate.markdown(joinpath(@__DIR__,"..","examples",example_source), output_dir,documenter=false)    
end

generated_examples=joinpath.("examples",readdir(output))


clean=true
if isinteractive()
    clean=false
end

makedocs(
    sitename="VoronoiFVM.jl",
    modules = [VoronoiFVM],
    clean = clean,
    doctest = false,
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
        "Examples" => generated_examples
    ]
)

if !isinteractive()
    deploydocs(repo = "github.com/j-fu/VoronoiFVM.jl.git")
end


