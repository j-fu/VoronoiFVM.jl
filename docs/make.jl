using Documenter, VoronoiFVM, Literate

#
# Replace SOURCE_URL marker with github url of source
#
function replace_source_url(input,source_url)
    lines_in = collect(eachline(IOBuffer(input)))
    lines_out=IOBuffer()
    for line in lines_in
        println(lines_out,replace(line,"SOURCE_URL" => source_url))
    end
    return String(take!(lines_out))
end




function make_all()
    #
    # Generate Markdown pages from examples
    #
    example_jl_dir = joinpath(@__DIR__,"..","examples")
    example_md_dir  = joinpath(@__DIR__,"src","examples")
    
    for example_source in readdir(example_jl_dir)
        base,ext=splitext(example_source)
        if ext==".jl"
            source_url="https://github.com/j-fu/VoronoiFVM.jl/raw/master/examples/"*example_source
            preprocess(buffer)=replace_source_url(buffer,source_url)
            Literate.markdown(joinpath(@__DIR__,"..","examples",example_source),
                              example_md_dir,
                              documenter=false,
                              info=false,
                              preprocess=preprocess)
        end
    end
    generated_examples=vcat(["runexamples.md"],joinpath.("examples",readdir(example_md_dir)))
    
    makedocs(
        sitename="VoronoiFVM.jl",
        modules = [VoronoiFVM],
        clean = false, 
        doctest = false,
        authors = "J. Fuhrmann",
        repo="https://github.com/j-fu/VoronoiFVM.jl",
        pages=[ 
            "Home"=>"index.md",
            "changes.md",
            "method.md",
            "API Documentation" => [
                "physics.md",
                "system.md",
                "diffeq.md",
                "allindex.md",
            ],
            "Examples" => generated_examples
        ]
    )
    
    rm(example_md_dir,recursive=true)
    
    if !isinteractive()
        deploydocs(repo = "github.com/j-fu/VoronoiFVM.jl.git")
    end
end

make_all()
