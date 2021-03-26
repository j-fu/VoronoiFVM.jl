##################################################################
"""
$(TYPEDEF)

Structure holding data for finite volume system solution.
Information on species distribution is kept in dense
matrices, and the solution array is of type Array{2}.


Unlike in the SparseSystem, the system matrix handles exactly those
degrees of freedom which correspond to unknowns, and dummy 
degrees of freedom where unknowns are not defined. Handling
of the sparse matrix structures for the bookeeping of the unknowns
has less overhead, but additional dummy equations are added
to the system matrix.

$(TYPEDFIELDS)
"""
mutable struct DenseSystem{Tv,Ti, Tm} <: AbstractSystem{Tv,Ti, Tm}
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
    Full matrix containing species numbers for inner regions
    """
    region_species::Array{Int8,2}

    """
    Full matrix containing species numbers for boundary regions
    """
    bregion_species::Array{Int8,2}

    """
    Full matrix containing degree of freedom numbers for each node
    """
    node_dof::Array{Int8,2}

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
    update::DenseSolutionArray{Tv}

    """
    Solution vector holding Newton residual
    """
    residual::DenseSolutionArray{Tv}

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

    DenseSystem{Tv,Ti, Tm}() where {Tv,Ti, Tm} = new()
end
##################################################################
"""
$(TYPEDSIGNATURES)

Constructor for DenseSystem.
"""
function  DenseSystem(grid,physics::Physics;matrixindextype=Int32)
    Tv=coord_type(grid)
    Ti=index_type(grid)
    Tm=matrixindextype
    this=DenseSystem{Tv,Ti,Tm}()
    maxspec=physics.num_species
    this.grid=grid
    this.physics=physics
    this.region_species=spzeros(Int8,Int16,maxspec,num_cellregions(grid))
    this.bregion_species=spzeros(Int8,Int16,maxspec,num_bfaceregions(grid))
    this.node_dof=spzeros(Int8,Int32,maxspec,num_nodes(grid))
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
num_dof(this::DenseSystem)= length(this.node_dof)



##################################################################
"""
$(TYPEDSIGNATURES)

Create a solution vector for dense system.
If inival is not specified, the entries of the returned vector are undefined.
"""
unknowns(sys::DenseSystem{Tv};inival=undef) where Tv = unknowns(Tv,sys,inival=inival)

##################################################################
"""
$(TYPEDSIGNATURES)

Create a solution vector for dense system.
If inival is not specified, the entries of the returned vector are undefined.
"""
function unknowns(Tu::Type,sys::DenseSystem{Tv};inival=undef) where Tv
    a=Array{Tu}(undef,size(sys.node_dof)...)
    if inival!=undef
        fill!(a,inival)
    end
    return a
end


Base.reshape(v::DenseSolutionArray{Tu},sys::DenseSystem{Tv}) where {Tu,Tv}=v


##################################################################
"""
$(SIGNATURES)

Reshape vector to fit as solution to system.
"""
function Base.reshape(v::AbstractVector{Tu}, sys::DenseSystem{Tv}) where {Tu,Tv}
    @assert  length(v)==num_dof(sys)
    nspec=num_species(sys)
    reshape(v,Int64(nspec),Int64(length(v)/nspec))
end


#
# Assemble dummy equations for inactive species
#
function _eval_and_assemble_inactive_species(this::DenseSystem,U,Uold,F)
    if this.species_homogeneous
        return
    end
    for inode=1:size(this.node_dof,2)
        for ispec=1:size(this.node_dof,1)
            if !isdof(this,ispec,inode)
                F[ispec,inode]+= U[ispec,inode]-Uold[ispec,inode];
                idof=dof(F,ispec,inode)
                this.matrix[idof,idof]+=1.0
            end
        end
    end
end

#
# Initialize values in inactive dof for dense system
#
function _initialize_inactive_dof!(U::DenseSolutionArray{Tv},this::DenseSystem{Tv}) where {Tv}
    if this.species_homogeneous
        return
    end
    for inode=1:size(this.node_dof,2)
        for ispec=1:size(this.node_dof,1)
            if !isdof(this,ispec,inode)
                U[ispec,inode]=0
            end
        end
    end
end


"""
   
"""
function LinearAlgebra.norm(sys::DenseSystem,u,p)
    _initialize_inactive_dof!(u,sys)
    norm(u,p)
end
