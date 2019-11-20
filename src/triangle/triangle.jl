module Triangle
using DocStringExtensions
#
# Prepare shared lib calls
#
const basedir=Base.@__DIR__
const depsdir=basedir*"/../../deps/"

# if Sys.iswindows()
# 	libsuffix = ".dll"

if Sys.isapple()
    libsuffix = ".dylib"
elseif Sys.islinux()
    libsuffix = ".so"
else
    Base.@error "Operating system not supported."
end

const libtriangle = depsdir*"usr/lib/libtriangle"*libsuffix


if ~isfile(libtriangle)
    Base.@error("Triangle library not found. Please run `Pkg.build(\"VoronoiFVM\")` first.")
end

#
# Struct maping Triangle's triangulateio
#
mutable struct CTriangulateIO
    pointlist :: Ptr{Cdouble}
    pointattributelist :: Ptr{Cdouble}
    pointmarkerlist :: Ptr{Cint}
    numberofpoints :: Cint
    numberofpointattributes :: Cint
    
    trianglelist :: Ptr{Cint}
    triangleattributelist :: Ptr{Cdouble}
    trianglearealist :: Ptr{Cdouble}
    neighborlist :: Ptr{Cint}
    numberoftriangles :: Cint
    numberofcorners :: Cint
    numberoftriangleattributes :: Cint
    
    segmentlist :: Ptr{Cint}
    segmentmarkerlist :: Ptr{Cint}
    numberofsegments :: Cint

    holelist :: Ptr{Cdouble}
    numberofholes :: Cint

    regionlist :: Ptr{Cdouble}
    numberofregions :: Cint

    edgelist :: Ptr{Cint}
    edgemarkerlist :: Ptr{Cint}
    normlist :: Ptr{Cdouble}
    numberofedges :: Cint
end 

#
# Constructor for `CTriangulateIO`. Initialize everything as `NULL`. Only for internal use.
# 
function CTriangulateIO()
    return CTriangulateIO(C_NULL, C_NULL, C_NULL, 0, 0,
                          C_NULL, C_NULL, C_NULL, C_NULL, 0, 0, 0,
                          C_NULL, C_NULL, 0,
                          C_NULL, 0,
                          C_NULL, 0,
                          C_NULL, C_NULL, C_NULL, 0)
end


#
# "Raw" triangulation call to Triangle library
#
function triangulate(triangle_switches::String, ctio_in::CTriangulateIO, ctio_out::CTriangulateIO, vor_out::CTriangulateIO)
    ccall((:triangulate,libtriangle),
          Cvoid,
          ( Cstring,
            Ref{CTriangulateIO}, 
            Ref{CTriangulateIO},
            Ref{CTriangulateIO}
            ),
          triangle_switches,
          Ref(ctio_in),
          Ref(ctio_out),
          Ref(vor_out)
          )
end



##################################################################################
"""
$(TYPEDEF)

Julia version of Triangle's triangulateio structure.

$(TYPEDFIELDS)
"""
mutable struct TriangulateIO
    pointlist :: Array{Cdouble,2}
    pointattributelist :: Array{Cdouble,2}
    
    trianglelist :: Array{Cint,2}
    triangleattributelist :: Array{Cdouble,2}
    trianglearealist :: Array{Cdouble,1} # input only
    neighborlist :: Array{Cint,2} # output only
    numberofcorners :: Cint
    
    segmentlist :: Array{Cint,2}
    segmentmarkerlist :: Array{Cint,1}

    holelist :: Array{Cdouble,2} # input only, copied to output list

    regionlist ::  Array{Cdouble,2} # input only, copied to output list

    edgelist :: Array{Cint,2} # output only
    edgemarkerlist :: Array{Cint,1} # output only
    normlist :: Array{Cdouble,2} # output only
end 

##########################################################
"""
$(TYPEDSIGNATURES)


Create empty TriangulateIO structure
"""
function TriangulateIO()
    return TriangulateIO(Array{Cdouble,2}(undef,0,0), # poinlist
                         Array{Cdouble,2}(undef,0,0), # pointattrlist
                         Array{Cint,2}(undef,0,0),   # trianglelist
                         Array{Cint,2}(undef,0,0),   # triangleattrlist
                         Array{Cdouble,1}(undef,0),   # trianglearealist
                         Array{Cint,2}(undef,0,0),   # nblist
                         0,
                         Array{Cint,2}(undef,0,0),   # seglist
                         Array{Cint,1}(undef,0),     # segmarkers
                         Array{Cdouble,2}(undef,0,0), # holelist
                         Array{Cdouble,2}(undef,0,0), # regionlist
                         Array{Cint,2}(undef,0,0),   # edgelist
                         Array{Cint,1}(undef,0),     # edgemarkerlist
                         Array{Cdouble,2}(undef,0,0)  # normlist
                         )
end


#
# Construct CTriangulateIO from TriangulateIO
#
function CTriangulateIO(tio::TriangulateIO)
    ctio=CTriangulateIO()
    ctio.numberofpoints=size(tio.pointlist,2)
    @assert ctio.numberofpoints>0
    @assert  size(tio.pointlist,1)==2
    ctio.pointlist=pointer(tio.pointlist)
    ctio.numberofpointattributes=size(tio.pointattributelist,1)
    if ctio.numberofpointattributes>0
        @assert size(tio.pointattributelist,2)==ctio.numberofpoints
        ctio.pointattributelist=pointer(tio.pointattributelist)
    end
    
    ctio.numberoftriangles=size(tio.trianglelist,2)
    if ctio.numberoftriangles>0
        ctio.trianglelist=pointer(tio.trianglelist)
        ctio.numberoftriangleattributes=size(tio.triangleattributelist,1)
        if ctio.numberoftriangleattributes>0
            @assert size(tio.triangleattributelist,2)==ctio.numberoftriangles
            ctio.triangleattributes=pointer(tio.triangleattributes)
        end
        if size(tio.trianglearealist,1)>0
            @assert size(tio.trianglearealist,1)==ctio.numberoftriangles
            ctio.trianglearealist=pointer(tio.trianglearealist)
        end
        ctio.numberofcorners=tio.numberofcorners
    end
    
    ctio.numberofsegments=size(tio.segmentlist,2)
    if ctio.numberofsegments>0
        ctio.segmentlist=pointer(tio.segmentlist)
        @assert size(tio.segmentmarkerlist,1)==ctio.numberofsegments
        ctio.segmentmarkerlist=pointer(tio.segmentmarkerlist)
    end
    
    ctio.numberofholes=size(tio.holelist,2)
    if ctio.numberofholes>0
        @assert size(tio.holelist,1)==2
        ctio.holeist=pointer(tio.holelist)
    end
    
    ctio.numberofregions=size(tio.regionlist,2)
    if ctio.numberofregions>0
        @assert size(tio.regionlist,1)==4
        ctio.regionlist=pointer(tio.regionlist)
    end
    return ctio
end


#
# Construct TriangulateIO from CTriangulateIO
#
function TriangulateIO(ctio::CTriangulateIO)
    tio=TriangulateIO()
    if ctio.numberofpoints>0
        tio.pointlist = convert(Array{Cdouble,2}, Base.unsafe_wrap(Array, ctio.pointlist, (2,Int(ctio.numberofpoints)), own=true))
    end
    if ctio.numberofpointattributes>0
        tio.pointattributelist=convert(Array{Cdouble,2}, Base.unsafe_wrap(Array, ctio.pointattributelistlist, (Int(ctio.numberofpointattributes),Int(ctio.numberofpoints)), own=true))
    end
    if ctio.numberoftriangles>0
        tio.trianglelist=convert(Array{Cint,2}, Base.unsafe_wrap(Array, ctio.trianglelist, (3,Int(ctio.numberoftriangles)), own=true))
    end
    if ctio.numberoftriangleattributes>0
        tio.triangleattributelist=convert(Array{Cdouble,2}, Base.unsafe_wrap(Array, ctio.triangleattributelist, (Int(ctio.numberoftriangleattributes),Int(ctio.numberoftriangles)), own=true))
    end
    # todo: trianglearealist check if the C pointer is 0
    if ctio.numberofsegments>0
        tio.segmentlist=convert(Array{Cint,2}, Base.unsafe_wrap(Array, ctio.segmentlist, (2,Int(ctio.numberofsegments)), own=true))
        tio.segmentmarkerlist=convert(Array{Cint,1}, Base.unsafe_wrap(Array, ctio.segmentmarkerlist, (Int(ctio.numberofsegments)), own=true))
    end
    
    # todo: copy regions,holes to output
    if ctio.numberofedges>0
        tio.edgelist=convert(Array{Cint,2}, Base.unsafe_wrap(Array, ctio.edgelist, (2,Int(ctio.numberofedges)), own=true))
        tio.edgemarkerlist=convert(Array{Cint,1}, Base.unsafe_wrap(Array, ctio.edgemarkerlist, (Int(ctio.numberofedges)), own=true))
        # todo: check if the C pointer is 0
        # tio.normlist=convert(Array{Cdouble,2}, Base.unsafe_wrap(Array, ctio.edgemarkerlist, (2,Int(ctio.numberofedges)), own=true))
    end
    return tio
end


##########################################################
"""
$(TYPEDSIGNATURES)

Create triangulation. Returns tuple `(tri_out::TriangulateIO, vor_out::TriangulateIO)`
containing the output triangulation and the optional Voronoi tesselation.
"""
function triangulate(triangle_switches::String, tio_in::TriangulateIO)
    ctio_in=CTriangulateIO(tio_in)
    ctio_out=CTriangulateIO()
    cvor_out=CTriangulateIO()
    triangulate(triangle_switches,ctio_in,ctio_out,cvor_out)
    tio_out=TriangulateIO(ctio_out)
    vor_out=TriangulateIO(cvor_out)
    return tio_out,vor_out
end


end # module Triangle
