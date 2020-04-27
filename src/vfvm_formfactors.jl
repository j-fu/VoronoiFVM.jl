#=
Form factor + edge velocity calculation
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



################################################
"""
$(SIGNATURES)

Calculate node volume  and voronoi surface contributions for cell.
""" 

function cellfactors!(::Type{Edge1D},::Type{Cartesian1D},coord,cellnodes,icell::Int,nodefac::Vector{Tv},edgefac::Vector{Tv}) where Tv
    K=cellnodes[1,icell]
    L=cellnodes[2,icell]
    xK=coord[1,K]
    xL=coord[1,L]
    d=abs(xL[1]-xK[1])
    nodefac[1]=d/2
    nodefac[2]=d/2
    edgefac[1]=1/d
    nothing
end

function cellfactors!(::Type{Edge1D},::Type{<:Polar1D}, coord,cellnodes,icell::Int,nodefac::Vector{Tv},edgefac::Vector{Tv}) where Tv
    K=cellnodes[1,icell]
    L=cellnodes[2,icell]
    xK=coord[1,K]
    xL=coord[1,L]
    r0=xK[1]
    r1=xL[1]
    if r1<r0
        r0=xL[1]
        r1=xK[1]
    end
    rhalf=0.5*(r1+r0);
    πv::Tv=π
    # cpar[1]= πv*(r1*r1-r0*r0);         # circular volume
    nodefac[1]= πv*(rhalf*rhalf-r0*r0);   # circular volume between midline and boundary
    nodefac[2]= πv*(r1*r1-rhalf*rhalf);   # circular volume between midline and boundary
    edgefac[1]= 2.0*πv*rhalf/(r1-r0);     # circular surface / width
    nothing
end
    
    
function cellfactors!(::Type{Triangle2D},::Type{<:Cartesian2D},coord,cellnodes,icell,npar,epar)
    i1=cellnodes[1,icell]
    i2=cellnodes[2,icell]
    i3=cellnodes[3,icell]

    # Fill matrix of edge vectors
    V11= coord[1,i2]- coord[1,i1]
    V21= coord[2,i2]- coord[2,i1]
    
    V12= coord[1,i3]- coord[1,i1]
    V22= coord[2,i3]- coord[2,i1]
    
    V13= coord[1,i3]- coord[1,i2]
    V23= coord[2,i3]- coord[2,i2]
    
    
    # Compute determinant 
    det=V11*V22 - V12*V21
    vol=0.5*det
    
    ivol = 1.0/vol
    
    # squares of edge lengths
    dd1=V13*V13+V23*V23 # l32
    dd2=V12*V12+V22*V22 # l31
    dd3=V11*V11+V21*V21 # l21
    
    
    # contributions to \sigma_kl/h_kl
    epar[1]= (dd2+dd3-dd1)*0.125*ivol
    epar[2]= (dd3+dd1-dd2)*0.125*ivol
    epar[3]= (dd1+dd2-dd3)*0.125*ivol
    
    
    # contributions to \omega_k
    npar[1]= (epar[3]*dd3+epar[2]*dd2)*0.25
    npar[2]= (epar[1]*dd1+epar[3]*dd3)*0.25
    npar[3]= (epar[2]*dd2+epar[1]*dd1)*0.25
    nothing
end                              


function cellfactors!(::Type{Triangle2D},::Type{<:Cylindrical2D},coord,cellnodes,icell::Int,npar::Vector{Tv},epar::Vector{Tv}) where Tv
    function area2d(coord1, coord2, coord3)
        V11= coord2[1]- coord1[1]
        V21= coord2[2]- coord1[2]
        
        V12= coord3[1]- coord1[1]
        V22= coord3[2]- coord1[2]
        
        V13= coord3[1]- coord2[1]
        V23= coord3[2]- coord2[2]
        
        # Compute determinant 
        det=V11*V22 - V12*V21
        area=abs(0.5*det)
    end
        

    πv::Tv=π
    i1=cellnodes[1,icell]
    i2=cellnodes[2,icell]
    i3=cellnodes[3,icell]
    
    
    # Fill matrix of edge vectors
    V11= coord[1,i2]- coord[1,i1]
    V21= coord[2,i2]- coord[2,i1]
    
    V12= coord[1,i3]- coord[1,i1]
    V22= coord[2,i3]- coord[2,i1]
    
    V13= coord[1,i3]- coord[1,i2]
    V23= coord[2,i3]- coord[2,i2]
    
    # Compute determinant 
    det=V11*V22 - V12*V21
    area=0.5*det
    
    # Integrate R over triangle (via quadrature rule)
    vol=2.0*πv*area*(coord[1,i1]+coord[1,i2]+coord[1,i3])/3.0
    
    # squares of edge lengths
    dd1=V13*V13+V23*V23 # l32
    dd2=V12*V12+V22*V22 # l31
    dd3=V11*V11+V21*V21 # l21
    
    emid23=[0.5*(coord[1,i3]+coord[1,i2]),
            0.5*(coord[2,i3]+coord[2,i2])]
    
    emid13=[0.5*(coord[1,i1]+coord[1,i3]),
            0.5*(coord[2,i1]+coord[2,i3])]
    
    emid12=[0.5*(coord[1,i1]+coord[1,i2]),
            0.5*(coord[2,i1]+coord[2,i2])]
    
    cc=Vector{Float64}(undef,2) # TODO: replace this allocation + views
    tricircumcenter!(cc,coord[:,i1],coord[:,i2],coord[:,i3])
    
    r(p)=p[1]
    z(p)=p[2]
    sq(x)=x*x
    epar[1]= πv*(r(cc)+r(emid23))*sqrt(sq(r(cc)-r(emid23))+sq(z(cc)-z(emid23)))/sqrt(dd1);
    epar[2]= πv*(r(cc)+r(emid13))*sqrt(sq(r(cc)-r(emid13))+sq(z(cc)-z(emid13)))/sqrt(dd2);
    epar[3]= πv*(r(cc)+r(emid12))*sqrt(sq(r(cc)-r(emid12))+sq(z(cc)-z(emid12)))/sqrt(dd3);
    
    rintegrate(coord1, coord2, coord3)=2.0*πv*area2d(coord1,coord2,coord3)*(coord1[1]+coord2[1]+coord3[1])/3.0
    npar[1]=rintegrate(coord[:,i1],cc,emid13)+rintegrate(coord[:,i1],cc,emid12)
    npar[2]=rintegrate(coord[:,i2],cc,emid23)+rintegrate(coord[:,i2],cc,emid12)
    npar[3]=rintegrate(coord[:,i3],cc,emid13)+rintegrate(coord[:,i3],cc,emid23)
    nothing
end                              


################################################
"""
$(SIGNATURES)

Calculate node volume  contributions for boundary face.
""" 

function bfacefactors!(::Type{Vertex0D},::Type{Cartesian1D},coord,bfacenodes,ibface::Int,nodefac::Vector{Tv}) where Tv
    nodefac[1]=1.0
    nothing
end

function bfacefactors!(::Type{Vertex0D},::Type{<:Polar1D},coord,bfacenodes,ibface::Int,nodefac::Vector{Tv}) where Tv
    inode=bfacenodes[1,ibface]
    r=coord[1,inode]
    nodefac[1]=2*pi*r
    nothing
end

# TODO: Test
function bfacefactors!(::Type{Edge1D},::Type{<:Cartesian2D},coord,bfacenodes,ibface::Int,nodefac::Vector{Tv}) where Tv
    i1=bfacenodes[1,ibface]
    i2=bfacenodes[2,ibface]
    dx=coord[1,i1]-coord[1,i2]
    dy=coord[2,i1]-coord[2,i2]
    d=0.5*sqrt(dx*dx+dy*dy)
    nodefac[1]=d
    nodefac[2]=d
    nothing
end

# TODO: Test
function bfacefactors!(::Type{Edge1D},::Type{<:Cylindrical2D},coord,bfacenodes,ibface::Int,nodefac::Vector{Tv}) where Tv
    i1=bfacenodes[1,ibface]
    i2=bfacenodes[2,ibface]
    r1=coord[1,i1]
    r2=coord[1,i2]
    z1=coord[2,i1]
    z2=coord[2,i2]
    dr=r1-r2
    rmid=(r1+r2)/2
    dz=z1-z2
    l=sqrt(dr*dr+dz*dz)
    nodefac[1]=pi*(r1+rmid)*l/2
    nodefac[2]=pi*(r2+rmid)*l/2
    nothing
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

    cn=grid[CellNodes]
    ec=grid[EdgeCells]
    en=grid[EdgeNodes]
    coord=grid[Coordinates]
    
    velovec=zeros(Float64,num_edges(grid))
    if dim_space(grid)==1
        for iedge=1:num_edges(grid)
            K=en[1,iedge]
            L=en[2,iedge]
            elen=coord[1,L]-coord[1,K]
            vx,vy=velofunc((coord[1,K]+coord[1,L])/2,0)
            velovec[iedge]=-elen*vx
        end
    else
        for iedge=1:num_edges(grid)
            K=en[1,iedge]
            L=en[2,iedge]
            p1=Vector{Float64}(undef,2)
            p2=Vector{Float64}(undef,2)
            tricircumcenter!(p1,
                             coord[:,cn[1,ec[1,iedge]]],
                             coord[:,cn[2,ec[1,iedge]]],
                             coord[:,cn[3,ec[1,iedge]]])
            if ec[2,iedge]>0
                tricircumcenter!(p2,
                                 coord[:,cn[1,ec[2,iedge]]],
                                 coord[:,cn[2,ec[2,iedge]]],
                                 coord[:,cn[3,ec[2,iedge]]])
            else
                p2.=0.5*(coord[:,K]+coord[:,L])
            end    
            hnormal=coord[:,K]-coord[:,L]
            velovec[iedge]=integrate(p1,p2,hnormal,velofunc)
        end
    end
    return velovec
end
