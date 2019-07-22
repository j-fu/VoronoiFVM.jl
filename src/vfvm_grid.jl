##########################################################
"""
$(TYPEDEF)

Abstract type for grid like datastructures [`VoronoiFVM.Grid`](@ref) and [`VoronoiFVM.SubGrid`](@ref).
"""
abstract type AbstractGrid end



##########################################################
"""
$(TYPEDSIGNATURES)

Space dimension of grid
"""
dim_space(grid::AbstractGrid)= size(grid.coord,1)


##########################################################
"""
$(TYPEDSIGNATURES)


Number of nodes in grid
"""
num_nodes(grid::AbstractGrid)= size(grid.coord,2)


##########################################################
"""
$(TYPEDSIGNATURES)

Number of cells in grid
"""
num_cells(grid::AbstractGrid)= size(grid.cellnodes,2)

##########################################################
"""
$(TYPEDSIGNATURES)

Return index of i-th local node in cell icell
"""
cellnode(grid::AbstractGrid,inode,icell)=grid.cellnodes[inode,icell]

##########################################################
"""
$(TYPEDSIGNATURES)

Return view of coordinates of node `inode`.
"""
nodecoord(grid::AbstractGrid,inode)=view(grid.coord,:,inode)

##########################################################
"""
$(TYPEDSIGNATURES)

Return number of nodes per cell in grid.
"""
num_nodes_per_cell(grid::AbstractGrid)= size(grid.cellnodes,1)

##########################################################
"""
$(TYPEDSIGNATURES)

Return element type of grid coordinates.
"""
Base.eltype(grid::AbstractGrid)=Base.eltype(grid.coord)

##########################################################
"""
$(TYPEDEF)

Structure holding grid data. It is parametrised by the
type Tc of coordinates.

$(TYPEDFIELDS)

"""
struct Grid{Tc} <: AbstractGrid

    """ 
    2D Array of node coordinates
    """
    coord::Array{Tc,2}

    
    """
    2D Array of node indices per grid cell
    """
    cellnodes::Array{Int32,2}

    
    """
    Array of cell region numbers
    """
    cellregions::Array{Int32,1}

    
    """
    2D Array of node indices per boundary face
    """
    bfacenodes::ElasticArray{Int32,2,1}      

    
    """
    Array of boundary face region numbers
    """
    bfaceregions::Array{Int32,1}

    
    """
    Number of inner cell regions. Stored in an array in order
    to keep the struct immutable.
    """
    num_cellregions::Array{Int32,1}

    
    """
    Number of boundary face  regions. Stored in an array in order
    to keep the struct immutable.
    """
    num_bfaceregions::Array{Int32,1}

    
    """
    2D Array describing local scheme of distributions nodes per cell edge.
    """
    celledgenodes::Array{Int32,2}

end

function Base.show(io::IO,grid::Grid{Tc}) where Tc
    str=@sprintf("%s(dim_space=%d, num_nodes=%d, num_cells=%d, num_bfaces=%d)",
                 typeof(grid),dim_space(grid),num_nodes(grid), num_cells(grid), num_bfaces(grid))
    println(io,str)
end

##########################################################
"""
$(SIGNATURES)

Constructor for 1D grid.

Construct 1D grid from an array of node cordinates.
It creates two boundary regions with index 1 at the left end and
index 2 at the right end.

Primal grid holding unknowns: marked by `o`, dual
grid marking control volumes: marked by `|`.

```@raw html
 o-----o-----o-----o-----o-----o-----o-----o-----o
 |--|-----|-----|-----|-----|-----|-----|-----|--|
```

"""
function Grid(X::Array{Tc,1}) where Tc

    coord=reshape(X,1,length(X))
    cellnodes=zeros(Int32,2,length(X)-1)
    cellregions=zeros(Int32,length(X)-1)
    for i=1:length(X)-1 
        cellnodes[1,i]=i
        cellnodes[2,i]=i+1
        cellregions[i]=1
    end
    bfacenodes=ElasticArray{Int32}(undef,1,2)
    bfaceregions=zeros(Int32,2)
    bfacenodes[1,1]=1
    bfacenodes[1,2]=length(X)
    bfaceregions[1]=1
    bfaceregions[2]=2
    num_cellregions=[maximum(cellregions)]
    num_bfaceregions=[maximum(bfaceregions)]
    celledgenodes=reshape(Int32[1 2],:,1)
    return Grid{Tc}(coord,
                    cellnodes,
                    cellregions,
                    bfacenodes,
                    bfaceregions,
                    num_cellregions,
                    num_bfaceregions,
                    celledgenodes)
end


##########################################################
"""
$(SIGNATURES)

Constructor for 2D grid
from coordinate arrays. 
Boundary region numbers count counterclockwise:

| location  |  number |
| --------- | ------- |
| south     |       1 |
| east      |       2 |
| north     |       3 |
| west      |       4 |

"""
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
                
        if (geq(x1,coord[1,n1],coord[1,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=4
            return
        end
        if (leq(xn,coord[1,n1],coord[1,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=2
            return
        end
        if (geq(y1,coord[2,n1],coord[2,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=1
            return
        end
        if (leq(yn,coord[2,n1],coord[2,n2]))
            ibface=ibface+1
            bfacenodes[1,ibface]=n1
            bfacenodes[2,ibface]=n2
	    bfaceregions[ibface]=3
            return
        end
    end
    
    
    num_nodes=nx*ny
    num_cells=2*(nx-1)*(ny-1)
    num_bfacenodes=2*(nx-1)+2*(ny-1)
    
    coord=zeros(Tc,2,num_nodes)
    cellnodes=zeros(Int32,3,num_cells)
    cellregions=zeros(Int32,num_cells)
    bfacenodes=ElasticArray{Int32}(undef,2,num_bfacenodes)
    resize!(bfacenodes,2,num_bfacenodes)
    bfaceregions=zeros(Int32,num_bfacenodes)
    
    ipoint=0
    for iy=1:ny
        for ix=1:nx
            ipoint=ipoint+1
            coord[1,ipoint]=X[ix]
            coord[2,ipoint]=Y[iy]
        end
    end
    @assert(ipoint==num_nodes)
    
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
    @assert(icell==num_cells)
    
    #lazy way to  create boundary grid

    ibface=0
    for icell=1:num_cells
        n1=cellnodes[1,icell]
	n2=cellnodes[2,icell]
	n3=cellnodes[3,icell]
        check_insert_bface(n1,n2)
	check_insert_bface(n1,n3)
	check_insert_bface(n2,n3)
    end
    @assert(ibface==num_bfacenodes)

    num_cellregions=[maximum(cellregions)]
    num_bfaceregions=[maximum(bfaceregions)]
    celledgenodes=[2 1 1 ;
                   3 3 2]

    return Grid{Tc}(coord,
                    cellnodes,
                    cellregions,
                    bfacenodes,
                    bfaceregions,
                    num_cellregions,
                    num_bfaceregions,
                    celledgenodes)
end




######################################################
"""
$(TYPEDSIGNATURES)

Edit region numbers of grid cells via rectangular mask.
"""
function cellmask!(grid::Grid,
                   maskmin::AbstractArray,
                   maskmax::AbstractArray,
                   ireg::Int;
                   eps=1.0e-10)
    xmaskmin=maskmin.-eps
    xmaskmax=maskmax.+eps
    for icell=1:num_cells(grid)
        in_region=true
        for inode=1:num_nodes_per_cell(grid)
            coord=nodecoord(grid,cellnode(grid,inode,icell))
            for idim=1:dim_space(grid)
                if coord[idim]<xmaskmin[idim]
                    in_region=false
                elseif coord[idim]>xmaskmax[idim]
                    in_region=false
                end
            end
        end
        if in_region
            grid.cellregions[icell]=ireg
        end
    end
    grid.num_cellregions[1]=max(num_cellregions(grid),ireg)
end


######################################################
"""
$(TYPEDSIGNATURES)

Edit region numbers of grid  boundary facets  via rectangular mask.
Currently, only for 1D grids, inner boundaries can be added.
"""
function bfacemask!(grid::Grid,
                    maskmin::AbstractArray,
                    maskmax::AbstractArray,
                    ireg::Int;
                    eps=1.0e-10)

    
    
    xmaskmin=maskmin.-eps
    xmaskmax=maskmax.+eps
    
    function isbface(ix)
        for ibface=1:num_bfaces(grid)
            if grid.bfacenodes[1,ibface]==ix
                return ibface
            end
            return 0
        end
    end
    if dim_space(grid)==1
        for inode=1:num_nodes(grid)
            x=grid.coord[1,inode]
            if x>xmaskmin[1] && x<xmaskmax[1]
                ibface=isbface(inode)
                if ibface>0
                    grid.bfaceregions[ibface]=ireg
                else
                    ibface=length(grid.bfaceregions)+1
                    push!(grid.bfaceregions,ireg)
                    append!(grid.bfacenodes,[inode])
                end
            end
        end
    else
        for ibface=1:num_bfaces(grid)
            in_region=true
            for inode=1:num_nodes_per_bface(grid)
                coord=nodecoord(grid,bfacenode(grid,inode,ibface))
                for idim=1:dim_space(grid)
                    if coord[idim]<xmaskmin[idim]
                        in_region=false
                    elseif coord[idim]>xmaskmax[idim]
                        in_region=false
                    end
                end
            end
            if in_region
                grid.bfaceregions[ibface]=ireg
            end
        end
    end
        
    grid.num_bfaceregions[1]=max(num_bfaceregions(grid),ireg)
end





################################################
"""
$(SIGNATURES)

Calculate node volume  and voronoi surface contributions for cell.
""" 
function cellfactors!(grid::Grid{Tv},icell::Int,nodefac::Vector{Tv},edgefac::Vector{Tv}) where Tv

    function cellfac1d!(grid::Grid{Tv},icell::Int,nodefac::Vector{Tv},edgefac::Vector{Tv}) where Tv
        K=cellnode(grid,1,icell)
        L=cellnode(grid,2,icell)
        xK=nodecoord(grid,K)
        xL=nodecoord(grid,L)
        d=abs(xL[1]-xK[1])
        nodefac[1]=d/2
        nodefac[2]=d/2
        edgefac[1]=1/d
    end
    
    function cellfac2d!(grid::Grid{Tv},icell::Int,npar::Vector{Tv},epar::Vector{Tv}) where Tv
        i1=cellnode(grid,1,icell)
        i2=cellnode(grid,2,icell)
        i3=cellnode(grid,3,icell)
        
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
    
    
    if dim_space(grid)==1
        cellfac1d!(grid,icell,nodefac,edgefac)
    elseif dim_space(grid)==2
        cellfac2d!(grid,icell,nodefac,edgefac)
    end
end

################################################
"""
$(SIGNATURES)

Calculate node volume  and voronoi surface contributions for boundary face.
""" 
function bfacefactors!(grid::Grid{Tv},icell::Int,nodefac::Vector{Tv}) where Tv

    # 1D bface form factors
    function bfacefac1d!(grid::Grid,ibface::Int,nodefac::Vector{Tv}) where Tv
        nodefac[1]=1.0
    end
    
    
    # 2D bface form factors
    function bfacefac2d!(grid::Grid,ibface::Int,nodefac::Vector{Tv}) where Tv
        i1=bfacenode(grid,1,ibface)
        i2=bfacenode(grid,2,ibface)
        dx=grid.coord[1,i1]-grid.coord[1,i2]
        dy=grid.coord[2,i1]-grid.coord[2,i2]
        d=0.5*sqrt(dx*dx+dy*dy)
        nodefac[1]=d
        nodefac[2]=d
    end
    
    
    if dim_space(grid)==1
        bfacefac1d!(grid,icell,nodefac)
    elseif dim_space(grid)==2
        bfacefac2d!(grid,icell,nodefac)
    end
end

################################################
"""
$(TYPEDSIGNATURES)

Bulk region number for cell
"""
reg_cell(grid::Grid,icell)=grid.cellregions[icell]

################################################
"""
$(TYPEDSIGNATURES)


Boundary region number for boundary face
"""
reg_bface(grid::Grid,icell)=grid.bfaceregions[icell]

################################################
"""
$(TYPEDSIGNATURES)


Topological dimension of grid
"""
dim_grid(grid::Grid)= size(grid.bfacenodes,1)

################################################
"""
$(TYPEDSIGNATURES)

Index of boundary face node.
"""
bfacenode(grid::Grid,inode,icell)=grid.bfacenodes[inode,icell]

################################################
"""
$(TYPEDSIGNATURES)

Index of cell edge node.
"""
celledgenode(grid::Grid,inode,iedge,icell)=grid.cellnodes[grid.celledgenodes[inode,iedge],icell]

################################################
"""
$(TYPEDSIGNATURES)
    
Number of edges per grid cell.
"""
num_edges_per_cell(grid::Grid)= size(grid.celledgenodes,2)

################################################
"""
$(TYPEDSIGNATURES)

Number of nodes per boundary face
"""
num_nodes_per_bface(grid::Grid)= size(grid.bfacenodes,1)

################################################
"""
$(TYPEDSIGNATURES)

Number of boundary faces in grid.
"""
num_bfaces(grid::Grid)= size(grid.bfacenodes,2)

################################################
"""
$(TYPEDSIGNATURES)

Number of cell regions in grid.
"""
num_cellregions(grid::Grid)= grid.num_cellregions[1]

################################################
"""
$(TYPEDSIGNATURES)

Number of boundary face regions in grid.
"""
num_bfaceregions(grid::Grid)=grid.num_bfaceregions[1]



##################################################################
"""
$(TYPEDEF)
    
Subgrid of parent grid (mainly for visualization purposes). Intended
to hold support of species which are not defined everywhere.

$(TYPEDFIELDS)
"""
struct SubGrid{Tc} <: AbstractGrid


    """
    Parent Grid
    """
    parent::Grid

    
    """
    Incidence between subgrid node numbers and node numbers
    in parent.
    """
    node_in_parent::Array{Int32,1}

    
    """ 
    2D Array of coordinates per grid node
    """
    coord::Array{Tc,2}

    
    """
    2D Array of node numbers per grid cell
    """
    cellnodes::Array{Int32,2}

    
end


##################################################################
# Default transform for subgrid creation
function _copytransform!(a::AbstractArray,b::AbstractArray)
    for i=1:length(a)
        a[i]=b[i]
    end
end

##################################################################
"""
$(TYPEDSIGNATURES)

Create subgrid of list of regions.
"""
function subgrid(parent::Grid,
                 subregions::AbstractArray;
                 transform::Function=_copytransform!,
                 boundary=false)
    Tc=Base.eltype(parent)
    
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
        sub_gdim=dim_grid(parent)-1
    else
        xregions=parent.cellregions
        xnodes=parent.cellnodes
        sub_gdim=dim_grid(parent)
    end
    
    nodemark=zeros(Int32,num_nodes(parent))
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
        transform(localcoord[:,inode],parent.coord[:,sub_nip[inode]])
    end
    
    return SubGrid(parent,sub_nip,localcoord,sub_cellnodes)
end


