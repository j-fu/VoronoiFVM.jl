using Documenter, VoronoiFVM, Literate

# Turn block comments into "normal" hash comments
# as they currently are not handled by Literate.jl.
function hashify_block_comments(input)
    lines_in = collect(eachline(IOBuffer(input)))
    lines_out=IOBuffer()
    line_number=0
    in_block_comment_region=false
    for line in lines_in
        line_number+=1
        if occursin(r"^#=", line)
            if in_block_comment_region
                error("line $(line_number): already in block comment region\n$(line)")
            end
            println(lines_out,replace(line,r"^#=" => "#"))
            in_block_comment_region=true
        elseif occursin(r"^=#", line)
            if !in_block_comment_region
                error("line $(line_number): not in block comment region\n$(line)")
            end
            println(lines_out,replace(line,r"^=#" => "#"))
            in_block_comment_region=false
        else
            if in_block_comment_region
                println(lines_out,"# "*line)
            else
                println(lines_out,line)
            end
        end
    end
    return String(take!(lines_out))
end


function make_all()
    #
    # Generate Markdown pages from examples
    #
    output_dir  = joinpath(@__DIR__,"src","examples")
    example_dir = joinpath(@__DIR__,"..","examples")
    
    for example_source in readdir(example_dir)
        base,ext=splitext(example_source)
        if ext==".jl"
            Literate.markdown(joinpath(@__DIR__,"..","examples",example_source),
                              output_dir,
                              documenter=false,
                              info=false,
                              preprocess=hashify_block_comments)
        end
    end
    generated_examples=joinpath.("examples",readdir(output_dir))
    
    makedocs(
        sitename="VoronoiFVM.jl",
        modules = [VoronoiFVM],
        clean = true,
        doctest = false,
        authors = "J. Fuhrmann",
        repo="https://github.com/j-fu/VoronoiFVM.jl",
        pages=[ 
            "Home"=>"index.md",
            "changes.md",
            "method.md",
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
end

make_all()
