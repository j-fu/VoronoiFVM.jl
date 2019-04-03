using SparseArrays

##################################################################
# Grid

abstract type AbstractGrid end


mutable struct GridData
    ncellregions::Int32
    nbfaceregions::Int32
    nedges_per_cell::Int32
    celledgenodes::Array{Int32,2}
    GridData()=new()
end

struct Grid{Tc} <: AbstractGrid
    nodecoord::Array{Tc,2} # point coordinates
    cellnodes::Array{Int32,2} 
    cellregions::Array{Int32,1}
    bfacenodes::Array{Int32,2} 
    bfaceregions::Array{Int32,1}
    griddata::GridData
    cellfactors::Function
    bfacefactors::Function
end




"""
    mutable struct Node

Structure holding local node information.
Fields:

    index::Int32
    coord::Array{Float64,1}

"""
mutable struct Node
    index::Int32
    coord::Array{Float64,1}
    region::Int32
    nspecies::Int32
    Node()=new()
end


"""
    mutable struct Edge

Structure holding local edge information.

Fields:

    index::Int32
    nodeK::Int32
    nodeL::Int32
    coordK::Array{Float64,1}
    coordL::Array{Float64,1}


"""
mutable struct Edge
    index::Int32
    nodeK::Int32
    nodeL::Int32
    region::Int32
    coordK::Array{Float64,1}
    coordL::Array{Float64,1}
    nspecies::Int32
    Edge()=new()
end

"""
   function edgelength(edge::Edge)
   
Calculate the length of an edge. 
"""
function edgelength(edge::Edge)
    l::Float64
    l=0.0
    for i=1:length(edge.coordK)
        d=edge.coordK[i]-edge.coordL[i]
        l=l+d*d
    end
    return l
end




function cellfac1d(grid,icell,nodefac,edgefac)
    K=cellnodes(grid,1,icell)
    L=cellnodes(grid,2,icell)
    xK=nodecoord(grid,K)
    xL=nodecoord(grid,L)
    d=abs(xL[1]-xK[1])
    nodefac[1]=d/2
    nodefac[2]=d/2
    edgefac[1]=1/d
end


function bfacefac1d(grid,ibface,nodefac)
    nodefac[1]=1.0
end

function cellfac2d(grid,icell,npar,epar)
    i1=cellnodes(grid,1,icell);
    i2=cellnodes(grid,2,icell)
    i3=cellnodes(grid,3,icell)
    p1=nodecoord(grid,i1)
    p2=nodecoord(grid,i2)
    p3=nodecoord(grid,i3)
    

    # Fill matrix of edge vectors
    V11= p2[1]- p1[1]
    V21= p2[2]- p1[2]
    
    V12= p3[1]- p1[1]
    V22= p3[2]- p1[2]
    
    V13= p3[1]- p2[1]
    V23= p3[2]- p2[2]
    
    
    
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
end                              


function bfacefac2d(grid,ibface,nodefac)
    i1=cellnodes(grid,1,ibface)
    i2=cellnodes(grid,2,ibface)
    p1=nodecoord(grid,i1)
    p2=nodecoord(grid,i2)
    dx=p1[1]-p2[1]
    dy=p1[2]-p2[2]
    d=0.5*sqrt(dx*dx+dy*dy)
    nodefac[1]=d
    nodefac[2]=d
end



function Grid(X::Array{Tc,1}) where Tc
    nodecoord=reshape(X,1,length(X))
    cellnodes=zeros(Int32,2,length(X)-1)
    cellregions=zeros(Int32,length(X)-1)
    for i=1:length(X)-1 
        cellnodes[1,i]=i
        cellnodes[2,i]=i+1
        cellregions[i]=1
    end
    bfacenodes=zeros(Int32,1,2)
    bfaceregions=zeros(Int32,2)
    bfacenodes[1,1]=1
    bfacenodes[1,2]=length(X)
    bfaceregions[1]=1
    bfaceregions[2]=2
    griddata=GridData()
    griddata.ncellregions=maximum(cellregions)
    griddata.nbfaceregions=maximum(bfaceregions)
    griddata.nedges_per_cell=1
    griddata.celledgenodes=reshape(Int32[1 2],:,1)
    return Grid{Tc}(nodecoord,
                    cellnodes,
                    cellregions,
                    bfacenodes,
                    bfaceregions,
                    griddata,
                    cellfac1d,bfacefac1d)
end


function  Grid(X::Array{Tc,1},Y::Array{Tc,1}) where Tc
    function leq(x, x1, x2)
        if (x>x1)
            return false
        end
        if (x>x2)
            return false
        end
        return true
    end
    
    function geq(x, x1, x2)
        if (x<x1)
            return false
        end
        if (x<x2)
            return false
        end
        return true
    end

    nx=length(X)
    ny=length(Y)
    
    hmin=X[2]-X[1]
    for i=1:nx-1
        h=X[i+1]-X[i]
        if h <hmin
            hmin=h
        end
    end
    for i=1:ny-1
        h=Y[i+1]-Y[i]
        if h <hmin
            hmin=h
        end
    end
    
    @assert(hmin>0.0)
    eps=1.0e-5*hmin

    x1=X[1]+eps
    xn=X[nx]-eps
    y1=Y[1]+eps
    yn=Y[ny]-eps
    
    
    function  check_insert_bface(n1,n2)
                
        if (geq(x1,nodecoord[1,n1],nodecoord[1,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=4
            return
        end
        if (leq(xn,nodecoord[1,n1],nodecoord[1,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=2
            return
        end
        if (geq(y1,nodecoord[2,n1],nodecoord[2,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=1
            return
        end
        if (leq(yn,nodecoord[2,n1],nodecoord[2,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=3
            return
        end
    end
    
    
    nnodes=nx*ny
    ncells=2*(nx-1)*(ny-1)
    nbfacenodes=2*(nx-1)+2*(ny-1)
    
    nodecoord=zeros(Tc,2,nnodes)
    cellnodes=zeros(Int32,3,ncells)
    cellregions=zeros(Int32,ncells)
    bfacenodes=zeros(Int32,2,nbfacenodes)
    bfaceregions=zeros(Int32,nbfacenodes)
    
    ipoint=0
    for iy=1:ny
        for ix=1:nx
            ipoint=ipoint+1
            nodecoord[1,ipoint]=X[ix]
            nodecoord[2,ipoint]=Y[iy]
        end
    end
    @assert(ipoint==nnodes)
    
    icell=0
    for iy=1:ny-1
        for ix=1:nx-1
	    ip=ix+(iy-1)*nx
	    p00 = ip
	    p10 = ip+1
	    p01 = ip  +nx
	    p11 = ip+1+nx
            
            icell=icell+1
            cellnodes[1,icell]=p00
            cellnodes[2,icell]=p10
            cellnodes[3,icell]=p11
            cellregions[icell]=1
            
            
            icell=icell+1
            cellnodes[1,icell]=p11
            cellnodes[2,icell]=p01
            cellnodes[3,icell]=p00
            cellregions[icell]=1
        end
    end
    @assert(icell==ncells)
    
    #lazy way to  create boundary grid

    ibface=0
    for icell=1:ncells
        n1=cellnodes[1,icell]
	n2=cellnodes[2,icell]
	n3=cellnodes[3,icell]
        check_insert_bface(n1,n2)
	check_insert_bface(n1,n3)
	check_insert_bface(n2,n3)
    end
    @assert(ibface==nbfacenodes)
    griddata=GridData()
    griddata.ncellregions=maximum(cellregions)
    griddata.nbfaceregions=maximum(bfaceregions)
    griddata.nedges_per_cell=3
    griddata.celledgenodes=[2 1 1 ;
                            3 3 2]
    return Grid{Tc}(nodecoord,cellnodes,
                    cellregions,bfacenodes,
                    bfaceregions,griddata,
                    cellfac2d,bfacefac2d)
end

function cellmask!(grid::Grid, maskmin::AbstractArray,maskmax::AbstractArray,ireg::Integer;eps=1.0e-10)
    xmaskmin=maskmin.-eps
    xmaskmax=maskmax.-eps
    for icell=1:ncells(grid)
        in_region=true
        for inode=1:nnodes_per_cell(grid)
            coord=nodecoord(grid,cellnodes(grid,inode,icell))
            for idim=1:spacedim(grid)
                if coord[idim]<maskmin[idim]
                    in_region=false
                elseif coord[idim]>maskmax[idim]
                    in_region=false
                end
            end
        end
        if in_region
            grid.cellregions[icell]=ireg
        end
    end
    grid.griddata.ncellregions=max(grid.griddata.ncellregions,ireg)
end



spacedim(grid::AbstractGrid)= size(grid.nodecoord,1)
nnodes(grid::AbstractGrid)= size(grid.nodecoord,2)
ncells(grid::AbstractGrid)= size(grid.cellnodes,2)
cellnodes(grid::AbstractGrid,inode,icell)=grid.cellnodes[inode,icell]
nodecoord(grid::AbstractGrid,inode)=view(grid.nodecoord,:,inode)
nnodes_per_cell(grid::AbstractGrid)= size(grid.cellnodes,1)


cellfactors(grid::Grid,icell,nodefac,edgefac)=grid.cellfactors(grid,icell,nodefac,edgefac)

bfacefactors(grid::Grid,icell,nodefac)=grid.bfacefactors(grid,icell,nodefac)

cellregions(grid::Grid,icell)=grid.cellregions[icell]

bfaceregions(grid::Grid,icell)=grid.bfaceregions[icell]

griddim(grid::Grid)= size(grid.bfacenodes,1)

bfacenodes(grid::Grid,inode,icell)=grid.bfacenodes[inode,icell]

celledgenodes(grid::Grid,inode,iedge,icell)=grid.cellnodes[grid.griddata.celledgenodes[inode,iedge],icell]

nedges_per_cell(grid::Grid)= grid.griddata.nedges_per_cell

nnodes_per_bface(grid::Grid)= size(grid.bfacenodes,1)

nbfaces(grid::Grid)= size(grid.bfacenodes,2)

ncellregions(grid::Grid)= grid.griddata.ncellregions

nbfaceregions(grid::Grid)=grid.griddata.nbfaceregions



##################################################################
# SubGrid
struct SubGrid{Tc} <: AbstractGrid
    parent::Grid
    cellnodes::Array{Int32,2}
    nodecoord::Array{Tc,2}
    node_in_parent::Array{Int32,1}
end


function copytransform!(a::AbstractArray,b::AbstractArray)
    for i=1:length(a)
        a[i]=b[i]
    end
end

function SubGrid(parent::Grid,subregions::AbstractArray;transform::Function=copytransform!,boundary=false)
    Tc=eltype(parent.nodecoord)
    
    @inline function insubregions(xreg)
        for i in eachindex(subregions)
            if subregions[i]==xreg
                return true
            end
        end
        return false
    end

    
    if boundary
        xregions=parent.bfaceregions
        xnodes=parent.bfacenodes
        sub_gdim=griddim(parent)-1
    else
        xregions=parent.cellregions
        xnodes=parent.cellnodes
        sub_gdim=griddim(parent)
    end
    
    nodemark=zeros(Int32,nnodes(parent))
    ncn=size(xnodes,1)
    
    nsubcells=0
    nsubnodes=0
    for icell in eachindex(xregions)
        if insubregions(xregions[icell])
            nsubcells+=1
            for inode=1:ncn
                ipnode=xnodes[inode,icell]
                if nodemark[ipnode]==0
                    nsubnodes+=1
                    nodemark[ipnode]=nsubnodes
                end
            end
        end
    end
    
    sub_cellnodes=zeros(Int32,ncn,nsubcells)
    sub_nip=zeros(Int32,nsubnodes)
    for inode in eachindex(nodemark)
        if nodemark[inode]>0
            sub_nip[nodemark[inode]]=inode
        end
    end
    
    isubcell=0
    for icell in eachindex(xregions)
        if insubregions(xregions[icell])
            isubcell+=1
            for inode=1:ncn
                ipnode=xnodes[inode,icell]
                sub_cellnodes[inode,isubcell]=nodemark[ipnode]
            end
        end
    end

    localcoord=zeros(Tc,sub_gdim,nsubnodes)
    @views for inode=1:nsubnodes
        transform(localcoord[:,inode],parent.nodecoord[:,sub_nip[inode]])
    end
    
    return SubGrid(parent,sub_cellnodes,localcoord,sub_nip)
end



