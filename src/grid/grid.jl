#=
Definition of grid structure and main constructors
=#


##########################################################
"""
    $(TYPEDEF)

Abstract type for grid like datastructures [`VoronoiFVM.Grid`](@ref) and [`VoronoiFVM.SubGrid`](@ref).
"""
abstract type AbstractGrid end





##########################################################

##########################################################
"""
$(TYPEDEF)

Type for dispatching formfactor calculations.
Intepret 1D coordinates as cartesian.
"""
struct Cartesian1D end

"""
$(TYPEDEF)

Type for dispatching formfactor calculations.
Intepret 2D coordinates as cartesian.
"""
struct Cartesian2D end

"""
$(TYPEDEF)

Type for dispatching formfactor calculations.
Intepret 3D coordinates as cartesian.
"""
struct Cartesian3D end

"""
$(TYPEDEF)

Type for dispatching formfactor calculations.
Intepret 1D coordinates as radial coordinate  r
assuming circular symmetry around origin 0.
"""
struct CircularSymmetric1D end

"""
$(TYPEDEF)

Type for dispatching formfactor calculations.
Intepret 2D coordinates as radial coordinate r
and eight coordinate z assuming circular symmetry 
around axis r=0.
"""
struct CircularSymmetric2D end

"""
$(TYPEDEF)

Type for dispatching formfactor calculations.
Intepret 1D coordinates as radial coordinate  r
assuming spehrical symmetry around point r=0.
"""
struct SphericalSymmetric1D end


##########################################################
"""
$(TYPEDEF)

Structure holding grid data. It is parametrised by the
type Tc of coordinates.

$(TYPEDFIELDS)

"""
mutable struct Grid{Tc,Ti} <: AbstractGrid

    """ 
    2D Array of node coordinates
    """
    coord::Array{Tc,2}

    
    """
    2D Array of node indices per grid cell
    """
    cellnodes::Array{Ti,2}

    """
    Array of edge indices per grid cell.
    Instatiated with prepare_edges!(grid)
    """
    celledges::Array{Ti,2}

    """
    Array of cell indices per grid edge.
    The number of these indices may be less than 
    the number of columns in this array. Nonexisting 
    indices are set to 0
    Instatiated with prepare_edges!(grid)
    """
    edgecells::Array{Ti,2}

    """
    Array of node indices per grid edge.
    Instatiated with prepare_edges!(grid)
    """
    edgenodes::Array{Ti,2}
    
    """
    Array of cell region numbers
    """
    cellregions::Array{Ti,1}

    
    """
    2D Array of node indices per boundary face
    """
    bfacenodes::Array{Ti,2}      

    
    """
    Array of boundary face region numbers
    """
    bfaceregions::Array{Ti,1}

    
    """
    Number of inner cell regions.
    """
    num_cellregions::Ti

    
    """
    Number of boundary face  regions.
    """
    num_bfaceregions::Ti

    
    """
    2D Array describing local scheme of distributions nodes per cell edge.
    """
    local_celledgenodes::Array{Ti,2}

    """
    Type of coordinate system
    """
    coord_type::DataType
end

function Base.show(io::IO,grid::Grid{Tc}) where Tc
    if num_edges(grid)>0
        str=@sprintf("%s(dim_space=%d, num_nodes=%d, num_cells=%d, num_bfaces=%d, num_edges=%d)",
                     typeof(grid),dim_space(grid),num_nodes(grid), num_cells(grid), num_bfaces(grid), num_edges(grid))
    else
        str=@sprintf("%s(dim_space=%d, num_nodes=%d, num_cells=%d, num_bfaces=%d)",
                     typeof(grid),dim_space(grid),num_nodes(grid), num_cells(grid), num_bfaces(grid))
        end
    println(io,str)
end


##########################################################
"""
  $(SIGNATURES)
        
  Main constructor for general grid from basic incidence information
"""
function Grid(coord::Array{Tc,2},
              cellnodes::Array{Ti,2},
              cellregions::Array{Ti,1},
              bfacenodes::Array{Ti,2},
              bfaceregions::Array{Ti,1}
              ) where {Tc,Ti}
    
    dim::Ti=size(coord,1)
    num_cellregions::Ti=maximum(cellregions)
    num_bfaceregions::Ti=maximum(bfaceregions)
    if dim==1
        local_celledgenodes=reshape(Ti[1 2],:,1)
        coord_type=Cartesian1D
    else
        # see grid/simplex.h in pdelib
        local_celledgenodes=zeros(Ti,2,3)
        local_celledgenodes[1,1]=2
        local_celledgenodes[2,1]=3

        local_celledgenodes[1,2]=3
        local_celledgenodes[2,2]=1
        
        local_celledgenodes[1,3]=1
        local_celledgenodes[2,3]=2

        coord_type=Cartesian2D
    end
        
    return Grid(coord,
                cellnodes,
                zeros(Ti,0,0),
                zeros(Ti,0,0),
                zeros(Ti,0,0),
                cellregions,
                bfacenodes,
                bfaceregions,
                num_cellregions,
                num_bfaceregions,
                local_celledgenodes,
                coord_type)
end

"""
  $(SIGNATURES)
        
  Set cartesian coordinates for grid.
"""
function cartesian!(grid)
    if dim_space(grid)==1
        grid.coord_type=Cartesian1D
    elseif dim_space(grid)==2
        grid.coord_type=Cartesian2D
    else dim_space(grid)==3
        grid.coord_type=Cartesian3D
    end
    return grid
end

"""
  $(SIGNATURES)
        
  Set circular coordinates for grid (1D or 2D).
"""
function circular_symmetric!(grid)
    if dim_space(grid)==1
        grid.coord_type=CircularSymmetric1D
    elseif dim_space(grid)==2
        grid.coord_type=CircularSymmetric2D
    else
        throw(DomainError(3,"Unable to handle circular symmetry for 3D grid"))
    end
    return grid
end

"""
  $(SIGNATURES)
        
  Set spherical coordinates for grid (1D).
"""
function spherical_symmetric!(grid)
    d=dim_space(grid)
    if d==1
        grid.coord_type=Spherical1D
    else
        throw(DomainError(d,"Unable to handle spherical symmetry for $(d)D grid"))
    end
    return grid
end



################################################
"""
$(SIGNATURES)

Prepare edge adjacencies (celledges, edgecells, edgenodes)
""" 
function prepare_edges!(grid)
    Ti=eltype(grid.cellnodes)
    
    # Create cell-node incidence matrix
    ext_cellnode_adj=ExtendableSparseMatrix{Ti,Ti}(num_nodes(grid),num_cells(grid))
    for icell=1:num_cells(grid)
        for inode=1:VoronoiFVM.num_nodes_per_cell(grid)
            ext_cellnode_adj[grid.cellnodes[inode,icell],icell]=1
        end
    end
    flush!(ext_cellnode_adj)
    # Get SparseMatrixCSC from the ExtendableMatrix
    cellnode_adj=ext_cellnode_adj.cscmatrix
    
    # Create node-node incidence matrix for neigboring
    # nodes. 
    nodenode_adj=cellnode_adj*transpose(cellnode_adj)

    # To get unique edges, we set the lower triangular part
    # including the diagonal to 0
    for icol=1:length(nodenode_adj.colptr)-1
        for irow=nodenode_adj.colptr[icol]:nodenode_adj.colptr[icol+1]-1
            if nodenode_adj.rowval[irow]>=icol
                nodenode_adj.nzval[irow]=0
            end
        end
    end
    dropzeros!(nodenode_adj)


    # Now we know the number of edges and
    nedges=length(nodenode_adj.nzval)

    
    if dim_space(grid)==2
        # Let us do the Euler test (assuming no holes in the domain)
        v=num_nodes(grid)
        e=nedges
        f=num_cells(grid)+1
        @assert v-e+f==2
    end
    if dim_space(grid)==1
        @assert nedges==num_cells(grid)
    end
    
    # Calculate edge nodes and celledges
    edgenodes=zeros(Ti,2,nedges)
    celledges=zeros(Ti,3,num_cells(grid))
    for icell=1:num_cells(grid)
        for iedge=1:VoronoiFVM.num_edges_per_cell(grid)
            n1=VoronoiFVM.celledgenode(grid,1,iedge,icell)
            n2=VoronoiFVM.celledgenode(grid,2,iedge,icell)            

            # We need to look in nodenod_adj for upper triangular part entries
            # therefore, we need to swap accordingly before looking
	    if (n1<n2)
		n0=n1
		n1=n2
		n2=n0;
	    end
            
            for irow=nodenode_adj.colptr[n1]:nodenode_adj.colptr[n1+1]-1
                if nodenode_adj.rowval[irow]==n2
                    # If the coresponding entry has been found, set its
                    # value. Note that this introduces a different edge orientation
                    # compared to the one found locally from cell data
                    celledges[iedge,icell]=irow
                    edgenodes[1,irow]=n1
                    edgenodes[2,irow]=n2
                end
            end
        end
    end


    # Create sparse incidence matrix for the cell-edge adjacency
    ext_celledge_adj=ExtendableSparseMatrix{Ti,Ti}(nedges,num_cells(grid))
    for icell=1:num_cells(grid)
        for iedge=1:VoronoiFVM.num_edges_per_cell(grid)
            ext_celledge_adj[celledges[iedge,icell],icell]=1
        end
    end
    flush!(ext_celledge_adj)
    celledge_adj=ext_celledge_adj.cscmatrix

    # The edge cell matrix is the transpose
    edgecell_adj=SparseMatrixCSC(transpose(celledge_adj))

    # Get the adjaency array from the matrix
    edgecells=zeros(Ti,2,nedges)
    for icol=1:length(edgecell_adj.colptr)-1
        ii=1
        for irow=edgecell_adj.colptr[icol]:edgecell_adj.colptr[icol+1]-1
            edgecells[ii,icol]=edgecell_adj.rowval[irow]
            ii+=1
        end
    end
    grid.edgecells=edgecells
    grid.celledges=celledges
    grid.edgenodes=edgenodes

    return grid
end






