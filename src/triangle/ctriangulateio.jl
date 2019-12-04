# Methods and structs from this file are not meant to be part ot the public interface, so
# documentation is in the comments, no docstrings here.
#

###############################################################
#
# Prepare shared library calls
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
    Base.@error("Triangle library not found. Please run `Pkg.build(\"TriangleRaw\")` first.")
end


###############################################################
# Check for LC_NUMERIC
# Possibly there are other options to handle this, e.g. replace strtod call in triangle
# See e.g.  https://github.com/JuliaLang/julia/issues/5928
# It might be the case that this problem occurs only in vendor distributions of Julia


function checklocale()
    endptr::Ptr{Ptr{UInt8}}=C_NULL
    if ccall((:strtod), Cdouble, (Cstring,Ptr{Ptr{UInt8}},), "0.5",endptr)!=0.5
        print(
            """ 

             Triangle uses  the C  library function strtod  to convert
             strings to double numbers.   The interpretation of double
             numbers  (whether  they  are  written with  ','  or  '.')
             depends on  the language  settings of your  computer.  In
             order  to ensure  portability  of  programs written  with
             TriangleRaw.jl, one should insist  on assuming '.' as the
             decimal point.

             In  the moment,  due to  whatever reason  (e.g. due  to a
             change  of  settings  performed  by  PyPlot),  the  wrong
             behavior has been detected.  The correct behaviour can be
             enforced by setting the environment variable 'LC_NUMERIC'
             to 'C' before starting Julia.

             See e.g. https://github.com/JuliaComputing/TextParse.jl/issues/30
         """)
        error("Missing or wrong value of LC_NUMERIC")
        return false
    end
    true
end

###############################################################
# Handling of triunsuitable callback

# Trivial default trinunsuitable function
function trivial_triunsuitable(org_x, org_y, dest_x, dest_y, apex_x, apex_y, area)
    return 0
end

# Global variable containing triunsuitable function
triunsuitable_func=trivial_triunsuitable

# triunsuitable function called from C (by triangulate(::ctriangulateio)) if -u flag has been set
function jl_wrap_triunsuitable(org_x::Cdouble, org_y::Cdouble, dest_x::Cdouble, dest_y::Cdouble, apex_x::Cdouble, apex_y::Cdouble, area::Cdouble)::Cint
    return Cint(triunsuitable_func(org_x, org_y, dest_x, dest_y, apex_x, apex_y, area))
end


#
# Struct mapping Triangle's triangulateio
#
# We use a mutable struct because we need fields to be able  initialized
# in differerent ways.
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
# Constructor for `CTriangulateIO`. Initializes everything as `NULL`. 
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
function triangulate(triangle_switches::String,
                     ctio_in::CTriangulateIO,
                     ctio_out::CTriangulateIO,
                     vor_out::CTriangulateIO)

    # Check locale settings for decimal point
    checklocale()

    # Set unsuitable callback
    if occursin("u",triangle_switches)
        c_wrap_triunsuitable=@cfunction(jl_wrap_triunsuitable, Cint, (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,))
        ccall((:triunsuitable_callback,libtriangle),Cvoid,(Ptr{Cvoid},),c_wrap_triunsuitable)
    end
    
    # Call triangulate
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


