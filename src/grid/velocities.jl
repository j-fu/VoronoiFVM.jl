#=
Handle edge velocity projections
=#

##################################################################
"""
$(SIGNATURES)

Find the circumcenter of a triangle.                 
									 
Created from C source of Jonathan R Shewchuk <jrs@cs.cmu.edu>

Modified to return absolute coordinates.
"""
function tricircumcenter!(circumcenter::Vector{Tv},a::Vector{Tv},b::Vector{Tv},c::Vector{Tv}) where Tv

    # Use coordinates relative to point `a' of the triangle.
    xba = b[1] - a[1]
    yba = b[2] - a[2]
    xca = c[1] - a[1]
    yca = c[2] - a[2]

    # Squares of lengths of the edges incident to `a'.
    balength = xba * xba + yba * yba
    calength = xca * xca + yca * yca
    
    # Calculate the denominator of the formulae.
    # if EXACT
    #    Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html	 
    #    to ensure a correctly signed (and reasonably accurate) result, 
    #    avoiding any possibility of division by zero.		    
    #  denominator = 0.5 / orient2d((double*) b, (double*) c, (double*) a)
    
    
    # Take your chances with floating-point roundoff
    denominator = 0.5 / (xba * yca - yba * xca)

    # Calculate offset (from `a') of circumcenter. 
    xcirca = (yca * balength - yba * calength) * denominator  
    ycirca = (xba * calength - xca * balength) * denominator  


# The result is returned both in terms of x-y coordinates and xi-eta	 
# coordinates, relative to the triangle's point `a' (that is, `a' is	 
# the origin of both coordinate systems).	 Hence, the x-y coordinates	 
# returned are NOT absolute; one must add the coordinates of `a' to	 
# find the absolute coordinates of the circumcircle.  However, this means	 
# that the result is frequently more accurate than would be possible if	 
# absolute coordinates were returned, due to limited floating-point	 
# precision.  In general, the circumradius can be computed much more	 
# accurately.								 
    

    circumcenter[1] = xcirca+a[1]
    circumcenter[2] = ycirca+a[2]

    return circumcenter
end




##################################################################
"""
$(SIGNATURES)

Project velocity onto grid edges,
"""
function edgevelocities(grid,velofunc)
    #
    # TODO: this should be generalized for more quadrules
    #
    function integrate(coordl,coordr,hnormal,velofunc)
        wl=1.0/6.0
        wm=2.0/3.0
        wr=1.0/6.0
        coordm=0.5*(coordl+coordr)
        (vxl,vyl)=velofunc(coordl[1],coordl[2])
        (vxm,vym)=velofunc(coordm[1],coordm[2])
        (vxr,vyr)=velofunc(coordr[1],coordr[2])
        return (wl*vxl +wm*vxm+wr*vxr)*hnormal[1] + (wl*vyl +wm*vym+wr*vyr)*hnormal[2]
    end
    
    @assert num_edges(grid)>0
    @assert dim_space(grid)<3

    cn=grid.cellnodes
    ec=grid.edgecells
    en=grid.edgenodes

    velovec=zeros(Float64,num_edges(grid))
    if VoronoiFVM.dim_space(grid)==1
        for iedge=1:num_edges(grid)
            K=en[1,iedge]
            L=en[2,iedge]
            elen=grid.coord[1,L]-grid.coord[1,K]
            vx,vy=velofunc((grid.coord[1,K]+grid.coord[1,L])/2,0)
            velovec[iedge]=-elen*vx
        end
    else
        for iedge=1:num_edges(grid)
            K=en[1,iedge]
            L=en[2,iedge]
            p1=Vector{Float64}(undef,2)
            p2=Vector{Float64}(undef,2)
            tricircumcenter!(p1,
                             grid.coord[:,cn[1,ec[1,iedge]]],
                             grid.coord[:,cn[2,ec[1,iedge]]],
                             grid.coord[:,cn[3,ec[1,iedge]]])
            if ec[2,iedge]>0
                tricircumcenter!(p2,
                                 grid.coord[:,cn[1,ec[2,iedge]]],
                                 grid.coord[:,cn[2,ec[2,iedge]]],
                                 grid.coord[:,cn[3,ec[2,iedge]]])
            else
                p2.=0.5*(grid.coord[:,K]+grid.coord[:,L])
            end    
            hnormal=grid.coord[:,K]-grid.coord[:,L]
            velovec[iedge]=integrate(p1,p2,hnormal,velofunc)
        end
    end
    return velovec
end
