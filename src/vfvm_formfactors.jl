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
function tricircumcenter!(circumcenter,a,b,c)

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

function cellfactors!(::Type{Edge1D},::Type{Cartesian1D},coord,cellnodes,icell,nodefac,edgefac)
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

function cellfactors!(::Type{Edge1D},::Type{<:Polar1D}, coord,cellnodes,icell,nodefac,edgefac)
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
    # cpar[1]= π*(r1*r1-r0*r0);         # circular volume
    nodefac[1]= π*(rhalf*rhalf-r0*r0);   # circular volume between midline and boundary
    nodefac[2]= π*(r1*r1-rhalf*rhalf);   # circular volume between midline and boundary
    edgefac[1]= 2.0*π*rhalf/(r1-r0);     # circular surface / width
    nothing
end
    
function cellfactors!(T::Type{Triangle2D},::Type{Cartesian2D},coord,cellnodes,icell,npar,epar)
#    cen=local_celledgenodes(T)
    
    @views n=SVector{3}(cellnodes[:,icell])
    V=@MMatrix zeros(2,3)
    dd=@MVector zeros(3)
    
    # Fill matrix of edge vectors
    V[1,3]= coord[1,n[2]]- coord[1,n[1]]
    V[2,3]= coord[2,n[2]]- coord[2,n[1]]
    
    V[1,2]= coord[1,n[3]]- coord[1,n[1]]
    V[2,2]= coord[2,n[3]]- coord[2,n[1]]
    
    V[1,1]= coord[1,n[3]]- coord[1,n[2]]
    V[2,1]= coord[2,n[3]]- coord[2,n[2]]
    
    
    # Compute determinant 
    det=V[1,3]*V[2,2] - V[1,2]*V[2,3]
    vol=0.5*det
    
    ivol = 1.0/vol
    
    # squares of edge lengths
    dd[3]=V[1,3]*V[1,3]+V[2,3]*V[2,3] # l32
    dd[2]=V[1,2]*V[1,2]+V[2,2]*V[2,2] # l31
    dd[1]=V[1,1]*V[1,1]+V[2,1]*V[2,1] # l21
    
    
    # contributions to \sigma_kl/h_kl
    epar[1]= (dd[2]+dd[3]-dd[1])*0.125*ivol
    epar[2]= (dd[3]+dd[1]-dd[2])*0.125*ivol
    epar[3]= (dd[1]+dd[2]-dd[3])*0.125*ivol
    
    
    # contributions to \omega_k
    npar[1]= (epar[3]*dd[3]+epar[2]*dd[2])*0.25
    npar[2]= (epar[1]*dd[1]+epar[3]*dd[3])*0.25
    npar[3]= (epar[2]*dd[2]+epar[1]*dd[1])*0.25

    nothing
end                              


function cellfactors!(::Type{Tetrahedron3D},::Type{Cartesian3D},coord,cellnodes,icell,npar,epar)
    # Transferred from WIAS/pdelib, (c) J. Fuhrmann, H. Langmach, I. Schmelzer
    es=@MVector zeros(6)
    dd=@MVector zeros(6)
    i1=cellnodes[1,icell]
    i2=cellnodes[2,icell]
    i3=cellnodes[3,icell]
    i4=cellnodes[4,icell]

    #-> celledgenodes
    pp1=(1,1,1,2,2,3)
    pp2=(2,3,4,3,4,4)
    
    pi1=(5,6,5,1)
    pi2=(6,3,1,4)
    pi3=(4,2,3,2)

    po1=(2,1,2,6)
    po2=(1,4,6,3)
    po3=(3,5,4,5)
    pedge=(6,5,4,3,2,1)

    for i=1:6
        p1=cellnodes[pp1[i],icell]
        p2=cellnodes[pp2[i],icell]
        dx=coord[1,p1]-coord[1,p2]
        dy=coord[2,p1]-coord[2,p2]
        dz=coord[3,p1]-coord[3,p2]
        dd[i]=dx*dx+dy*dy+dz*dz
        es[i]=0.0
        epar[i]=0.0
    end
    
    x1=coord[1,i2]-coord[1,i1]
    y1=coord[2,i2]-coord[2,i1]
    z1=coord[3,i2]-coord[3,i1]
    
    x2=coord[1,i3]-coord[1,i1]
    y2=coord[2,i3]-coord[2,i1]
    z2=coord[3,i3]-coord[3,i1]
    
    x3=coord[1,i4]-coord[1,i1]
    y3=coord[2,i4]-coord[2,i1]
    z3=coord[3,i4]-coord[3,i1]
    
    det= (x1*(y2*z3-y3*z2) + x2*(y3*z1-y1*z3) + x3*(y1*z2-y2*z1))

    if det <0
        det=-det
    end
    vol=det/6
    vv= 96*6*vol
    
    for  i=1:4
        npar[i]=0.0
        i1=pi1[i]
        i2=pi2[i]
        i3=pi3[i]
        
        h1=dd[i1]*(dd[i2]+dd[i3]-dd[i1])
        h2=dd[i2]*(dd[i3]+dd[i1]-dd[i2])
        h3=dd[i3]*(dd[i1]+dd[i2]-dd[i3])

        df=h1+h2+h3
        vf=(h1*dd[po1[i]] + h2*dd[po2[i]] + h3*dd[po3[i]]-2*dd[i1]*dd[i2]*dd[i3]) / (vv*df)
        es[i1] += h1*vf
        es[i2] += h2*vf
        es[i3] += h3*vf
    end
    
    for i=1:6
        epar[pedge[i]] = 6*es[i]/dd[i]
        npar[pp1[i]] += es[i]
        npar[pp2[i]] += es[i]
    end
    nothing
end


function cellfactors!(::Type{Triangle2D},::Type{<:Cylindrical2D},coord,cellnodes,icell,npar,epar)
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
    vol=2.0*π*area*(coord[1,i1]+coord[1,i2]+coord[1,i3])/3.0
    
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
    
    cc=@MVector zeros(2)
    tricircumcenter!(cc,coord[:,i1],coord[:,i2],coord[:,i3])
    
    r(p)=p[1]
    z(p)=p[2]
    sq(x)=x*x
    epar[1]= π*(r(cc)+r(emid23))*sqrt(sq(r(cc)-r(emid23))+sq(z(cc)-z(emid23)))/sqrt(dd1);
    epar[2]= π*(r(cc)+r(emid13))*sqrt(sq(r(cc)-r(emid13))+sq(z(cc)-z(emid13)))/sqrt(dd2);
    epar[3]= π*(r(cc)+r(emid12))*sqrt(sq(r(cc)-r(emid12))+sq(z(cc)-z(emid12)))/sqrt(dd3);
    
    rintegrate(coord1, coord2, coord3)=2.0*π*area2d(coord1,coord2,coord3)*(coord1[1]+coord2[1]+coord3[1])/3.0
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

function bfacefactors!(::Type{Vertex0D},::Type{Cartesian1D},coord,bfacenodes,ibface,nodefac,edgefac)
    nodefac[1]=1.0
    nothing
end

function bfacefactors!(::Type{Vertex0D},::Type{<:Polar1D},coord,bfacenodes,ibface,nodefac, edgefac)
    inode=bfacenodes[1,ibface]
    r=coord[1,inode]
    nodefac[1]=2*π*r
    nothing
end

function bfacefactors!(::Type{Edge1D},::Type{<:Cartesian2D},coord,bfacenodes,ibface,nodefac, edgefac)
    i1=bfacenodes[1,ibface]
    i2=bfacenodes[2,ibface]
    dx=coord[1,i1]-coord[1,i2]
    dy=coord[2,i1]-coord[2,i2]
    d=0.5*sqrt(dx*dx+dy*dy)
    nodefac[1]=d
    nodefac[2]=d
    nothing
end

function bfacefactors!(::Type{Edge1D},::Type{<:Cylindrical2D},coord,bfacenodes,ibface,nodefac,edgefac)
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
    nodefac[1]=π*(r1+rmid)*l/2
    nodefac[2]=π*(r2+rmid)*l/2
    nothing
end


function bfacefactors!(::Type{Triangle2D}, ::Type{<:Cartesian3D}, coord, bfacenodes, ibface, npar, epar)

    # Transferred from WIAS/pdelib, (c) J. Fuhrmann, H. Langmach, I. Schmelzer
    i1 = bfacenodes[1,ibface]
    i2 = bfacenodes[2,ibface]
    i3 = bfacenodes[3,ibface]
    dd = @MVector zeros(3)
    for j=1:3
        d = coord[j,i1] - coord[j,i3];  dd[2]+= d*d;
        d = coord[j,i2] - coord[j,i1];  dd[3]+= d*d;
        d = coord[j,i3] - coord[j,i2];  dd[1]+= d*d;
    end
    
    
    # Kanten-Flaechenanteile (ohne Abschneiden); epar als Hilfsfeld benutzt
    epar[1] = (dd[2]+dd[3]-dd[1])*dd[1];
    epar[2] = (dd[3]+dd[1]-dd[2])*dd[2];
    epar[3] = (dd[1]+dd[2]-dd[3])*dd[3];
    vol     = sqrt(epar[1]+epar[2]+epar[3])*0.25;
    
    d = 1.0/(8*vol);
        
    # Knoten-Flaechenanteile (ohne Abschneiden)
    npar[1] = (epar[3]+epar[2])*d*0.25;
    npar[2] = (epar[1]+epar[3])*d*0.25;
    npar[3] = (epar[2]+epar[1])*d*0.25;
    
    # Kantengewichte 
    epar[1] = epar[1]*d/dd[1];
    epar[2] = epar[2]*d/dd[2];
    epar[3] = epar[3]*d/dd[3];
    nothing
end

##################################################################

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


"""
$(SIGNATURES)

Project velocity onto grid edges,
"""
function edgevelocities(grid,velofunc)
    
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
            p1=@MVector zeros(2)
            p2=@MVector zeros(2)
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


function bfacevelocities(grid,velofunc)
    @assert dim_space(grid)<3
    bfacenodes=grid[BFaceNodes]
    coord=grid[Coordinates]
    bfacecells=grid[BFaceCells]
    bfacenormals=grid[BFaceNormals]
    bfr=grid[BFaceRegions]
    velovec=zeros(Float64,2,num_bfaces(grid))
    if dim_space(grid)==1
        for ibface=1:num_bfaces(grid)
            vx,vy=velofunc(coord[1,bfacenodes[1,ibface]])
            velovec[ibface]=vx*bfacenormals[1,ibface]
        end
    else
        for ibface=1:num_bfaces(grid)
            p1=coord[:,bfacenodes[1,ibface]]
            p2=coord[:,bfacenodes[2,ibface]]
            pm=0.5*(p1+p2)
            velovec[1,ibface]=integrate(p1,pm,bfacenormals[:,ibface],velofunc)
            velovec[2,ibface]=integrate(pm,p2,bfacenormals[:,ibface],velofunc)
        end
    end            
    return velovec
end
