#=
Grid constructors from various file formats
=#

######################################################
"""
    $(TYPEDSIGNATURES)
    
    Read grid from file.
"""
function Grid(::Type{<:IOStream};file::String="test.sg",format="")
    (fbase,fext)=splitext(file)
    if format==""
        format=fext[2:end]
    end
    @assert format=="sg"
    
    tks=TokenStream(file)
    expecttoken(tks,"SimplexGrid")
    version=parse(Float64,gettoken(tks))
    version20=false;

    if (version==2.0)
        version20=true;
    elseif (version==2.1)
        version20=false;
    else
        error("Read grid: wrong format version: $(version)")
    end

    dim::Int32=0
    coord=Array{Float64,2}(undef,0,0)
    cells=Array{Int32,2}(undef,0,0)
    regions=Array{Int32,1}(undef,0)
    faces=Array{Int32,2}(undef,0,0)
    bregions=Array{Int32,1}(undef,0)
    while(true)
        if (trytoken(tks,"DIMENSION"))
            dim=parse(Int32,gettoken(tks));
        elseif (trytoken(tks,"NODES")) 
            nnodes=parse(Int32,gettoken(tks));
            embdim=parse(Int32,gettoken(tks));
            if(dim!=embdim)
                error("Dimension error (DIMENSION $(dim)) in section NODES")
            end
            coord=Array{Float64,2}(undef,dim,nnodes)
            for inode=1:nnodes
                for idim=1:embdim
                    coord[idim,inode]=parse(Float64,gettoken(tks))
                end
            end
        elseif (trytoken(tks,"CELLS"))
            ncells=parse(Int32,gettoken(tks));
            cells=Array{Int32,2}(undef,dim+1,ncells)
            regions=Array{Int32,1}(undef,ncells)
            for icell=1:ncells
                for inode=1:dim+1
                    cells[inode,icell]=parse(Int32,gettoken(tks));
                end
                regions[icell]=parse(Int32,gettoken(tks));
	        if version20
		    for j=1:dim+1
		        gettoken(tks);  # skip file format garbage
                    end
                end
            end
        elseif (trytoken(tks,"FACES"))
            nfaces=parse(Int32,gettoken(tks));
            faces=Array{Int32,2}(undef,dim,nfaces)
            bregions=Array{Int32,1}(undef,nfaces)
            for iface=1:nfaces
                for inode=1:dim
                    faces[inode,iface]=parse(Int32,gettoken(tks));
                end
                bregions[iface]=parse(Int32,gettoken(tks));
	        if (version20)
		    for j=1:dim+2
		        gettoken(tks); #skip file format garbage
                    end
                end
            end
        else
            expecttoken(tks,"END")
            break
        end
    end
    Grid(coord,cells,regions,faces,bregions);
end
