using Documenter, VoronoiFVM, Literate, PlutoSliderServer, ExtendableGrids

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




function make_all(;with_examples=true)
    
    generated_examples=[]
    notebookmd=[]
    example_jl_dir = joinpath(@__DIR__,"..","examples")
    example_md_dir  = joinpath(@__DIR__,"src","examples")
    notebook_html_dir  = joinpath(@__DIR__,"src","nbhtml")

    ENV["RUNNING_DOCUMENTER"]=1
    if with_examples
        #
        # Run notebooks
        #
        notebooks=["api-update.jl",
                   "flux-reconstruction.jl",
                   "problemcase.jl"
                   ]
        
        # Use sliderserver to generate html
        
        export_directory(joinpath(@__DIR__,"..","pluto-examples"),
                          notebook_paths=notebooks,
                          Export_output_dir=joinpath(notebook_html_dir),
                          Export_offer_binder=false)

        # generate frame markdown for each notebook
        for notebook in notebooks
            base=split(notebook,".")[1]
            mdstring=
"""
# $(notebook)

Please note that in the html page, interactive elements like sliders are disabled.
Download this notebook: [$(notebook)](https://github.com/j-fu/VoronoiFVM.jl/blob/master/pluto-examples/$(notebook)).

```@raw html
<iframe style="height:15000px" width="100%" src="../$(base).html"> </iframe>
```
"""
            mdname=base*".md"
            push!(notebookmd,joinpath("nbhtml",mdname))
            io=open(joinpath(notebook_html_dir,mdname),"w")
            write(io,mdstring)
            close(io)
        end     
        #
        # Generate Markdown pages from examples
        #
        
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
    end

    
    makedocs(
        sitename="VoronoiFVM.jl",
        modules = [VoronoiFVM],
        clean = false, 
        doctest = true,
        authors = "J. Fuhrmann",
        repo="https://github.com/j-fu/VoronoiFVM.jl",
        pages=[ 
            "Home"=>"index.md",
            "changes.md",
            "method.md",
            "notebooks"=> notebookmd,
            "API Documentation" => [                
                "system.md",
                "physics.md",
                "solutions.md",
                "solver.md",
                "post.md",
                "misc.md",
                "quantities.md",
                "allindex.md",
            ],
            "Examples" => generated_examples
        ]
    )

    if with_examples
        rm(example_md_dir,recursive=true)
        rm(notebook_html_dir,recursive=true)
    end
    
    if !isinteractive()
        deploydocs(repo = "github.com/j-fu/VoronoiFVM.jl.git")
    end
end

make_all(with_examples=true)
