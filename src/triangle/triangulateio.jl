#
# Triangulateio structure and public interface methods
#
# Public interface methods have docstrings, internal ones
# comments.


#########################################################
"""
$(TYPEDEF)

Julia version of Triangle's triangulateio structure.


Used to  pass data into and  out of the triangulate()  procedure.  The
arrays are stored  in column major order with the first index being
the local number in a triangle or segment and the second index being the
index in the element count. This exactly corresponds
to  the layout  used in  Triangle. This  means that  all data  in this
structure are passed from/to triangle without copying.
                                                                         
Arrays are used to store points, triangles, markers, and so forth.       
                                                                         
Description of fields:

$(TYPEDFIELDS)

"""
mutable struct TriangulateIO

    """
    An array of point coordinates with `size(pointlist,1)==2`.
    `pointlist' must always point to a list of points. __Mandatory__.                   
    """
    pointlist :: Array{Cdouble,2}

    """
    An array of point attributes. There can be several attributes per point.
    Optional for input.
    """
    pointattributelist :: Array{Cdouble,2}

    """
    An array of point markers. Optional for input.
    """
    pointmarkerlist :: Array{Cint,1}

    """
    An array of triangle corners. The first three entries
    of each column describe the three nodes of the triangle
    in counterclockwise manner. They are followd by  any other nodes if the triangle    
    represents a nonlinear element.

    __Mandatory if the 'r' switch is used__. In this case `trianglelist' 
    must point to a list oftriangles with optional higher order nodes. 
    """    
    trianglelist :: Array{Cint,2}

    """
    An array of triangle attributes. There can be several attributes per triangle.
    Optional on input.
    """
    triangleattributelist :: Array{Cdouble,2}

    """
    An array of triangle area constraints. Input only.

    __Mandatory if both the 'r' and the 'a'  switch (with no number following) are used.__
    """
    trianglearealist :: Array{Cdouble,1}

    """
    An array of triangle neighbors. `size(neighborlist,1)==3`
    triangle.  Output only.
    """
    neighborlist :: Array{Cint,2}

    """
    An array of segment endpoints. `size(segmentlist,1)==2`

    __Mandatory if the 'p' switch is used__.
    """
    segmentlist :: Array{Cint,2}

    """
    An array of segment markers. Optional on input. If not set then
    segment markers on output default to zero.
    """
    segmentmarkerlist :: Array{Cint,1}

    """
    An array of holes. Holes are marked by some point from within the hole.
    Input only, although the array pointer is copied to the     
    output structure for convenience.

    Used if  the 'p' switch is used without the 'r' switch.
    """
    holelist :: Array{Cdouble,2}

    """
    An array of regional attributes and area constraints.     

    Used if  the 'p' switch is used without the 'r' switch.

    Each of the columns of this array contains 
    the constraint's x and y coordinates are at indices [1] and [2], 
    followed by the regional attribute at index [3], followed by the       
    maximum area at index [4]. So we have `size(regionlist,1)==4`.
    Note that each regional attribute is  
    used only if you select the 'A' switch, and each area constraint is    
    used only if you select the 'a' switch (with no number following), but 
    omitting one of these switches does not change the memory layout.      
    Input only, although the pointer is copied to the output structure for 
    convenience.                                                      
    """
    regionlist ::  Array{Cdouble,2}

    """
    An array of edge endpoints. `sizeof(edgelist,1)==2`.  Output only.
    """
    edgelist :: Array{Cint,2}

    """
    An array of edge markers; Output   only.
    """
    edgemarkerlist :: Array{Cint,1}

    """
    An array of normal vectors, used for infinite rays in       
    Voronoi diagrams. For eachfinite edge in a Voronoi diagram, 
    the normal vector written is the     
     zero vector.  `sizeof(normlist,1)==2`. Output only.                        
    """              
    normlist :: Array{Cdouble,2}
end 



##########################################################
"""
$(TYPEDSIGNATURES)

Return number of points in triangulatio structure.

"""
numberofpoints(tio::TriangulateIO)=size(tio.pointlist,2)


##########################################################
"""
$(TYPEDSIGNATURES)

Return number of segments in triangulateio structure.

"""
numberofsegments(tio::TriangulateIO)=size(tio.segmentlist,2)

##########################################################
"""
$(TYPEDSIGNATURES)

Return number of triangles in triangulateio structure.

"""
numberoftriangles(tio::TriangulateIO)=size(tio.trianglelist,2)

##########################################################
"""
$(TYPEDSIGNATURES)

Create TriangulateIO structure with empty data.

"""
function TriangulateIO()
    return TriangulateIO(Array{Cdouble,2}(undef,0,0), # poinlist
                         Array{Cdouble,2}(undef,0,0), # pointattrlist
                         Array{Cint,1}(undef,0),     # pointmarkers
                         Array{Cint,2}(undef,0,0),   # trianglelist
                         Array{Cint,2}(undef,0,0),   # triangleattrlist
                         Array{Cdouble,1}(undef,0),   # trianglearealist
                         Array{Cint,2}(undef,0,0),   # nblist
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

    if size(tio.pointmarkerlist,1)>0
        @assert size(tio.pointmarkerlist,1)==ctio.numberofpoints
        ctio.pointmarkerlist=pointer(tio.pointmarkerlist)
    end
    
    ctio.numberoftriangles=size(tio.trianglelist,2)
    if ctio.numberoftriangles>0
        ctio.trianglelist=pointer(tio.trianglelist)
        ctio.numberoftriangleattributes=size(tio.triangleattributelist,1)
        ctio.numberofcorners=size(tio.trianglelist,1) # JRS: "Each triangle occupies   `numberofcorners' ints."
        if ctio.numberoftriangleattributes>0
            @assert size(tio.triangleattributelist,2)==ctio.numberoftriangles
            ctio.triangleattributes=pointer(tio.triangleattributes)
        end
        if size(tio.trianglearealist,1)>0
            @assert size(tio.trianglearealist,1)==ctio.numberoftriangles
            ctio.trianglearealist=pointer(tio.trianglearealist)
        end
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

    if ctio.pointmarkerlist!=C_NULL
        tio.pointmarkerlist=convert(Array{Cint,1}, Base.unsafe_wrap(Array, ctio.pointmarkerlist, (Int(ctio.numberofpoints)), own=true))
    end


    if ctio.numberoftriangles>0
        tio.trianglelist=convert(Array{Cint,2}, Base.unsafe_wrap(Array, ctio.trianglelist, (3,Int(ctio.numberoftriangles)), own=true))
    end
    if ctio.numberoftriangleattributes>0
        tio.triangleattributelist=convert(Array{Cdouble,2}, Base.unsafe_wrap(Array, ctio.triangleattributelist, (Int(ctio.numberoftriangleattributes),Int(ctio.numberoftriangles)), own=true))
    end

    if ctio.trianglearealist!=C_NULL
        tio.trianglearealist=convert(Array{Cdouble,1}, Base.unsafe_wrap(Array, ctio.trianglearealist, (Int(ctio.numberoftriangles)), own=true))
    end

    if ctio.numberofsegments>0  && ctio.segmentlist!=C_NULL && ctio.segmentmarkerlist!=C_NULL
        tio.segmentlist=convert(Array{Cint,2}, Base.unsafe_wrap(Array, ctio.segmentlist, (2,Int(ctio.numberofsegments)), own=true))
        tio.segmentmarkerlist=convert(Array{Cint,1}, Base.unsafe_wrap(Array, ctio.segmentmarkerlist, (Int(ctio.numberofsegments)), own=true))
    end
    
    if ctio.numberofholes>0
        tio.holelist=convert(Array{Cdouble,2}, Base.unsafe_wrap(Array, ctio.holelist, (2,Int(ctio.numberofholes)), own=true))
    end
    
    if ctio.numberofregions>0
        tio.regionlist=convert(Array{Cdouble,2}, Base.unsafe_wrap(Array, ctio.holelist, (2,Int(ctio.numberofregions)), own=true))
    end
    
    if ctio.numberofedges>0
        tio.edgelist=convert(Array{Cint,2}, Base.unsafe_wrap(Array, ctio.edgelist, (2,Int(ctio.numberofedges)), own=true))
        tio.edgemarkerlist=convert(Array{Cint,1}, Base.unsafe_wrap(Array, ctio.edgemarkerlist, (Int(ctio.numberofedges)), own=true))
        if ctio.normlist!=C_NULL
            tio.normlist=convert(Array{Cdouble,2}, Base.unsafe_wrap(Array, ctio.normlist, (2,Int(ctio.numberofedges)), own=true))
        end
    end
    return tio
end


##########################################################
"""
$(TYPEDSIGNATURES)

Create triangulation. Returns tuple `(out::TriangulateIO, vor_out::TriangulateIO)`
containing the output triangulation and the optional Voronoi tesselation.

After a call to triangulate(), the valid fields of `out' and `vorout'
will depend, in an obvious way, on the choice of switches used.  Note
that when the 'p' switch is used, the pointers `holelist' and
`regionlist' are copied from `tr_in' to `out', but no new space is
#allocated;  On   the other hand, Triangle will never copy the `pointlist' pointer (or any 
others); new space is allocated for `out.pointlist', or if the `N'
switch is used, `out.pointlist' remains uninitialized. 

This is the list of switches used by triangle:

| Switch | Meaning                                                       |
|--------|:--------------------------------------------------------------|
| -p     | Triangulates a Planar Straight Line Graph (.poly file).       |
| -r     | Refines a previously generated mesh.                          |
| -q     | Quality mesh generation.  A minimum angle may be specified.   |
| -a     | Applies a maximum triangle area constraint.                   |
| -u     | Applies a user-defined triangle constraint.                   |
| -A     | Applies attributes to identify triangles in certain regions.  |
| -c     | Encloses the convex hull with segments.                       |
| -D     | Conforming Delaunay:  all triangles are truly Delaunay.       |
| -j     | Jettison unused vertices from output .node file.              |
| -e     | Generates an edge list.                                       |
| -v     | Generates a Voronoi diagram.                                  |
| -n     | Generates a list of triangle neighbors.                       |
| -g     | Generates an .off file for Geomview.                          |
| -B     | Suppresses output of boundary information.                    |
| -P     | Suppresses output of .poly file.                              |
| -N     | Suppresses output of .node file.                              |
| -E     | Suppresses output of .ele file.                               |
| -I     | Suppresses mesh iteration numbers.                            |
| -O     | Ignores holes in .poly file.                                  |
| -X     | Suppresses use of exact arithmetic.                           |
| -z     | Numbers all items starting from zero (rather than one).       |
| -o2    |  Generates second-order subparametric elements.               |
| -Y     | Suppresses boundary segment splitting.                        |
| -S     | Specifies maximum number of added Steiner points.             |
| -i     | Uses incremental method, rather than divide-and-conquer.      |
| -F     | Uses Fortune's sweepline algorithm, rather than d-and-c.      |
| -l     | Uses vertical cuts only, rather than alternating cuts.        |
| -s     | Force segments into mesh by splitting (instead of using CDT). |
| -C     | Check consistency of final mesh.                              |
| -Q     | Quiet:  No terminal output except errors.                     |
| -V     | Verbose:  Detailed information on what I'm doing.             |

"""
function triangulate(switches::String, tri_in::TriangulateIO)
    ctio_in=CTriangulateIO(tri_in)
    ctio_out=CTriangulateIO()
    cvor_out=CTriangulateIO()
    triangulate(switches,ctio_in,ctio_out,cvor_out)
    out=TriangulateIO(ctio_out)
    vor_out=TriangulateIO(cvor_out)
    return out,vor_out
end


##########################################################
"""
$(TYPEDSIGNATURES)

Set triunsuitable callback used by Triangle if the '-u' flag is set.

This is a function called by Triangle with the coordinates of the vertices
of a triangle in order to learn if that triangle needs further 
refinement (i.e. 'true' returned) or not ('false' returned).

Other checks (e.g. maximum edge lengths) are possible here as well.

Note, that the handling of this function is currently not thread safe.
````
function unsuitable(x1,y1,x2,y2,x3,y3,area)
    myarea=locally_desired_area(x1,y1,x2,y2,x3,y3)
    if area>myarea 
       return true
    else 
       return false
    end
end
````

"""
function triunsuitable(unsuitable::Function;check_signature=true)
    global triunsuitable_func
    if check_signature
        unsuitable(1,2,3,4,5,6,7)
    end
    triunsuitable_func=unsuitable
    nothing
end



