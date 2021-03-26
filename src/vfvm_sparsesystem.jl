##################################################################
"""
    $(TYPEDEF)
    
    Structure holding data for finite volume system solution.
Information on species distribution is kept in sparse
matrices, and the solution array is of type SparseSolutionArray,
i.e. effectively it is a sparse matrix.

Unlike in the DenseSystem, the system matrix handles exactly those
degrees of freedom which correspond to unknowns. However, handling
of the sparse matrix structures for the bookkeeping of the unknowns
creates overhead.

$(TYPEDFIELDS)
"""
mutable struct SparseSystem{Tv,Ti, Tm} <: AbstractSystem{Tv,Ti, Tm}

    """
    Grid
    """
    grid

    """
    Physics data
    """
    physics::Physics

    """
    Array of boundary values 
    """
    boundary_values::Array{Tv,2} 

    """
    Array of boundary factors 
    """
    boundary_factors::Array{Tv,2}

    """
    Sparse matrix containing species numbers for inner regions
    """
    region_species::SparseMatrixCSC{Int8,Ti}

    """
    Sparse matrix containing species numbers for boundary regions
    """
    bregion_species::SparseMatrixCSC{Int8,Ti}


    """
    Sparse matrix containing degree of freedom numbers for each node
    """
    node_dof::SparseMatrixCSC{Int8,Tm}

    """
    Jacobi matrix for nonlinear problem
    """
    matrix::ExtendableSparseMatrix{Tv,Tm}

    """
    Flag which says if the number of unknowns per node is constant
    """
    species_homogeneous::Bool

    """
    Solution vector holding Newton update
    """
    update::SparseSolutionArray{Tv,Tm}

    """
    Solution vector holding Newton residual
    """
    residual::SparseSolutionArray{Tv,Tm}

    """
    Precomputed geometry factors for cell nodes
    """
    cellnodefactors::Array{Tv,2}
    
    """
    Precomputed geometry factors for cell edges
    """
    celledgefactors::Array{Tv,2}

    """
    Precomputed geometry factors for boundary nodes
    """
    bfacenodefactors::Array{Tv,2}

    """
    Sparse matrix for generic operator handling
    """
    generic_matrix::SparseMatrixCSC

    """
    Sparse matrix colors for generic operator handling
    """
    generic_matrix_colors::Vector

    """
    Hash value of latest unknowns vector the assembly was called with
    """
    uhash::UInt64

    
    SparseSystem{Tv,Ti,Tm}() where {Tv,Ti,Tm} = new()
end

##################################################################

"""
$(TYPEDSIGNATURES)

Constructor for SparseSystem.
"""
function  SparseSystem(grid,physics::Physics; matrixindextype=Int32)
    Tv=coord_type(grid)
    Ti=index_type(grid)
    Tm=matrixindextype
    this=SparseSystem{Tv,Ti,Tm}()
    maxspec=physics.num_species
    this.grid=grid
    this.physics=physics
    this.region_species=spzeros(Int8,Int16,maxspec,num_cellregions(grid))
    this.bregion_species=spzeros(Int8,Int16,maxspec,num_bfaceregions(grid))
    this.node_dof=spzeros(Int8,Tm,maxspec,num_nodes(grid))
    this.boundary_values=zeros(Tv,maxspec,num_bfaceregions(grid))
    this.boundary_factors=zeros(Tv,maxspec,num_bfaceregions(grid))
    this.species_homogeneous=false
    this.uhash=0x0
    return this
end


##################################################################
"""
$(TYPEDSIGNATURES)

Number of degrees of freedom for system.
"""
num_dof(this::SparseSystem)= nnz(this.node_dof)



##################################################################
"""
$(SIGNATURES)

Create a solution vector for sparse system. 
The entries of the returned vector are undefined.
"""
unknowns(sys::SparseSystem{Tv,Ti,Tm};inival=undef) where {Tv,Ti, Tm}=unknowns(Tv,sys,inival=inival)

##################################################################
"""
$(SIGNATURES)

Create a solution vector for sparse system with given type.
If inival is not specified, the entries of the returned vector are undefined.
"""
function unknowns(Tu::Type, sys::SparseSystem{Tv, Ti, Tm};inival=undef) where {Tv,Ti,Tm}
    a0=Array{Tu}(undef,num_dof(sys))
    if inival!=undef
        fill!(a0,inival)
    end
    return SparseSolutionArray{Tu,Tm}(SparseMatrixCSC(sys.node_dof.m,
                                                      sys.node_dof.n,
                                                      sys.node_dof.colptr,
                                                      sys.node_dof.rowval,
                                                      a0
                                                      )
                                      )
end



##################################################################
"""
$(SIGNATURES)

Reshape vector to fit as solution to system.
"""
function Base.reshape(v::AbstractVector{Tu},sys::SparseSystem{Tv,Ti,Tm}) where {Tu,Tv,Ti,Tm}
    @assert  length(v)==num_dof(sys)
    SparseSolutionArray{Tu,Ti}(SparseMatrixCSC(sys.node_dof.m,
                                            sys.node_dof.n,
                                            sys.node_dof.colptr,
                                            sys.node_dof.rowval,
                                            Vector{Tu}(v)
                                            )
                            )
end

Base.reshape(v::SparseSolutionArray,sys::SparseSystem)=v




#
# Dummy routine for sparse system
#
function _eval_and_assemble_inactive_species(this::SparseSystem,U,Uold,F)
end

#
# Dummy routine for sparse system
#
function     _initialize_inactive_dof!(U::AbstractMatrix{Tv},this::SparseSystem{Tv}) where {Tv}
end

LinearAlgebra.norm(sys::SparseSystem,u,p)=norm(u.node_dof.nzval,p)

