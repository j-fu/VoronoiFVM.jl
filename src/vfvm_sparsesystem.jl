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
    Matrix factorization
    """
    factorization
    
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

    """
    Data for allocation check
    """
    allocs::Int
    
    SparseSystem{Tv,Ti,Tm}() where {Tv,Ti,Tm} = new()
end

##################################################################

"""
$(SIGNATURES)

Constructor for SparseSystem.
"""
function  SparseSystem(grid,physics::Physics; matrixindextype=Int32)
    Tv=coord_type(grid)
    Ti=index_type(grid)
    Tm=matrixindextype
    system=SparseSystem{Tv,Ti,Tm}()
    maxspec=physics.num_species
    system.grid=grid
    system.physics=physics
    system.region_species=spzeros(Int8,Int16,maxspec,num_cellregions(grid))
    system.bregion_species=spzeros(Int8,Int16,maxspec,num_bfaceregions(grid))
    system.node_dof=spzeros(Int8,Tm,maxspec,num_nodes(grid))
    system.boundary_values=zeros(Tv,maxspec,num_bfaceregions(grid))
    system.boundary_factors=zeros(Tv,maxspec,num_bfaceregions(grid))
    system.species_homogeneous=false
    system.uhash=0x0
    system.allocs=-1000
    system.factorization=nothing
    return system
end


##################################################################
"""
$(SIGNATURES)

Number of degrees of freedom for system.
"""
num_dof(system::SparseSystem)= nnz(system.node_dof)



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
function unknowns(Tu::Type, system::SparseSystem{Tv, Ti, Tm};inival=undef) where {Tv,Ti,Tm}
    a0=Array{Tu}(undef,num_dof(system))
    if inival!=undef
        fill!(a0,inival)
    end
    return SparseSolutionArray{Tu,Tm}(SparseMatrixCSC(system.node_dof.m,
                                                      system.node_dof.n,
                                                      system.node_dof.colptr,
                                                      system.node_dof.rowval,
                                                      a0
                                                      )
                                      )
end



##################################################################
"""
$(SIGNATURES)

Reshape vector to fit as solution to system.
"""
function Base.reshape(v::AbstractVector{Tu},system::SparseSystem{Tv,Ti,Tm}) where {Tu,Tv,Ti,Tm}
    @assert  length(v)==num_dof(system)
    SparseSolutionArray{Tu,Ti}(SparseMatrixCSC(system.node_dof.m,
                                            system.node_dof.n,
                                            system.node_dof.colptr,
                                            system.node_dof.rowval,
                                            Vector{Tu}(v)
                                            )
                            )
end

Base.reshape(v::SparseSolutionArray,sys::SparseSystem)=v




#
# Dummy routine for sparse system
#
function _eval_and_assemble_inactive_species(system::SparseSystem,U,Uold,F)
end

#
# Dummy routine for sparse system
#
function     _initialize_inactive_dof!(U::AbstractMatrix{Tv},system::SparseSystem{Tv}) where {Tv}
end

"""
$(SIGNATURES)

Calculate norm, paying attention to species distribution over regions
"""
LinearAlgebra.norm(system::SparseSystem,u,p)=norm(u.node_dof.nzval,p)

