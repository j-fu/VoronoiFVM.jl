#=
Form factor calculation
=#

################################################
"""
$(SIGNATURES)

Calculate node volume  and voronoi surface contributions for cell.
""" 
cellfactors!(grid::Grid{Tv},icell::Int,nodefac::Vector{Tv},edgefac::Vector{Tv}) where Tv=cellfactors!(grid.coord_type,grid, icell, nodefac,edgefac)


function cellfactors!(::Type{<:Cartesian1D},grid::Grid{Tv},icell::Int,nodefac::Vector{Tv},edgefac::Vector{Tv}) where Tv
    K::Int=cellnode(grid,1,icell)
    L::Int=cellnode(grid,2,icell)
    xK=nodecoord(grid,K)
    xL=nodecoord(grid,L)
    d=abs(xL[1]-xK[1])
    nodefac[1]=d/2
    nodefac[2]=d/2
    edgefac[1]=1/d
    nothing
end

function cellfactors!(::Type{<:CircularSymmetric1D}, grid::Grid{Tv},icell::Int,nodefac::Vector{Tv},edgefac::Vector{Tv}) where Tv
    K::Int=cellnode(grid,1,icell)
    L::Int=cellnode(grid,2,icell)
    xK=nodecoord(grid,K)
    xL=nodecoord(grid,L)
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
    
    
function cellfactors!(::Type{<:Cartesian2D},grid::Grid,icell,npar,epar)
    i1::Int=cellnode(grid,1,icell)
    i2::Int=cellnode(grid,2,icell)
    i3::Int=cellnode(grid,3,icell)
    coord=grid.coord

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


function cellfactors!(::Type{<:CircularSymmetric2D},grid::Grid{Tv},icell::Int,npar::Vector{Tv},epar::Vector{Tv}) where Tv
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
    i1::Int=cellnode(grid,1,icell)
    i2::Int=cellnode(grid,2,icell)
    i3::Int=cellnode(grid,3,icell)
    
    coord=grid.coord
    
    # Fill matrix of edge vectors
    V11= grid.coord[1,i2]- grid.coord[1,i1]
    V21= grid.coord[2,i2]- grid.coord[2,i1]
    
    V12= grid.coord[1,i3]- grid.coord[1,i1]
    V22= grid.coord[2,i3]- grid.coord[2,i1]
    
    V13= grid.coord[1,i3]- grid.coord[1,i2]
    V23= grid.coord[2,i3]- grid.coord[2,i2]
    
    # Compute determinant 
    det=V11*V22 - V12*V21
    area=0.5*det
    
    # Integrate R over triangle (via quadrature rule)
    vol=2.0*πv*area*(coord[1,i1]+coord[1,i2]+coord[1,i3])/3.0
    
    # squares of edge lengths
    dd1=V13*V13+V23*V23 # l32
    dd2=V12*V12+V22*V22 # l31
    dd3=V11*V11+V21*V21 # l21
    
    emid23=[0.5*(grid.coord[1,i3]+grid.coord[1,i2]),
            0.5*(grid.coord[2,i3]+grid.coord[2,i2])]
    
    emid13=[0.5*(grid.coord[1,i1]+grid.coord[1,i3]),
            0.5*(grid.coord[2,i1]+grid.coord[2,i3])]
    
    emid12=[0.5*(grid.coord[1,i1]+grid.coord[1,i2]),
            0.5*(grid.coord[2,i1]+grid.coord[2,i2])]
    
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
bfacefactors!(grid::Grid{Tv},icell::Int,nodefac::Vector{Tv}) where Tv=bfacefactors!(grid.coord_type,grid,icell,nodefac)

function bfacefactors!(::Type{<:Cartesian1D},grid::Grid,ibface::Int,nodefac::Vector{Tv}) where Tv
    nodefac[1]=1.0
    nothing
end

function bfacefactors!(::Type{<:CircularSymmetric1D},grid::Grid,ibface::Int,nodefac::Vector{Tv}) where Tv
    inode::Int=bfacenode(grid,1,ibface)
    r=grid.coord[1,i]
    nodefac[1]=2*pi*r
    nothing
end

# TODO: Test
function bfacefactors!(::Type{<:Cartesian2D},grid::Grid,ibface::Int,nodefac::Vector{Tv}) where Tv
    i1::Int=bfacenode(grid,1,ibface)
    i2::Int=bfacenode(grid,2,ibface)
    dx=grid.coord[1,i1]-grid.coord[1,i2]
    dy=grid.coord[2,i1]-grid.coord[2,i2]
    d=0.5*sqrt(dx*dx+dy*dy)
    nodefac[1]=d
    nodefac[2]=d
    nothing
end

# TODO: Test
function bfacefactors!(::Type{<:CircularSymmetric2D},grid::Grid,ibface::Int,nodefac::Vector{Tv}) where Tv
    i1::Int=bfacenode(grid,1,ibface)
    i2::Int=bfacenode(grid,2,ibface)
    r1=grid.coord[1,i1]
    r2=grid.coord[1,i2]
    z1=grid.coord[2,i1]
    z2=grid.coord[2,i2]
    dr=r1-r2
    rmid=(r1+r2)/2
    dz=z1-z2
    l=sqrt(dr*dr+dz*dz)
    nodefac[1]=pi*(r1+rmid)*l/2
    nodefac[2]=pi*(r2+rmid)*l/2
    nothing
end
