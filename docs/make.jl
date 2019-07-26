using Documenter, VoronoiFVM, Literate


#
# Generate Markdown pages from examples
#
output_dir  = joinpath(@__DIR__,"src","examples")
example_dir = joinpath(@__DIR__,"..","examples")
for example_source in readdir(example_dir)
    Literate.markdown(joinpath(@__DIR__,"..","examples",example_source), output_dir,documenter=false,info=false)    
end
generated_examples=joinpath.("examples",readdir(output_dir))

makedocs(
    sitename="VoronoiFVM.jl",
    modules = [VoronoiFVM],
    clean = !isinteractive(),
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


