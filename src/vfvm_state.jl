mutable struct SystemState{Tv, Ti, TSolArray, TData}
    system::VoronoiFVM.System

    """
    Parameter data 
    """
    data::TData

    """
    Solution vector
    """
    solution::TSolArray

    """
    Jacobi matrix for nonlinear problem
    """
    matrix::Union{ExtendableSparseMatrixCSC{Tv, Ti},
                  MTExtendableSparseMatrixCSC{Tv, Ti},
                  STExtendableSparseMatrixCSC{Tv, Ti},
                  Tridiagonal{Tv, Vector{Tv}},
                  #                  MultidiagonalMatrix,
                  BandedMatrix{Tv}}

    """
    Parameter derivative (vector of solution arrays)
    """
    dudp::Vector{TSolArray}
    
    """
    Vector holding Newton update
    """
    update::TSolArray

    """
    Vector holding Newton residual
    """
    residual::TSolArray

    """
    History record for last solution process
    """
    history::Any

    """
    Linear solver cache
    """
    linear_cache::Union{Nothing, LinearSolve.LinearCache}

    """
    Hash value of latest unknowns vector the assembly was called with
    """
    uhash::UInt64
end


"""
- `unknown_storage`: string or symbol.  
    Information  on  species  distribution  is kept  in  sparse  or  dense
    matrices matrices and, correspondingly, the  solution argray is of type
    SparseSolutionArray  or matrix,  respectively. In  the case  of sparse
    unknown storage,  the system matrix  handles exactly those  degrees of
    freedom which correspond to unknowns.  However, handling of the sparse
    matrix  structures  for  the   bookkeeping  of  the  unknowns  creates
    overhead.
     - `:dense` :  solution vector is an  `nspecies` x `nnodes`  dense matrix
     - `:sparse` :  solution vector is an `nspecies` x `nnodes`  sparse matrix
- `matrixindextype`: Integer type. Index type for sparse matrices created in the system.
"""
function SystemState(::Type{Tv}, system::AbstractSystem{Tv0, Tc, Ti, Tm}; data=system.physics.data) where {Tv,Tv0,Tc, Ti, Tm}
    lock(sysmutatelock)
    try
        _complete!(system)
        update_grid!(system)
    finally
        unlock(sysmutatelock)
    end

    nspec = size(system.node_dof, 1)
    n = num_dof(system)

    matrixtype = system.matrixtype
    #    matrixtype=:sparse
    # Sparse even in 1D is not bad, 

    if matrixtype == :default
        if !isdensesystem(system)
            matrixtype = :sparse
        else
            if nspec == 1
                matrixtype = :tridiagonal
            else
                matrixtype = :banded
            end
        end
    end

    if matrixtype == :tridiagonal
        matrix = Tridiagonal(zeros(Tv, n - 1), zeros(Tv, n), zeros(Tv, n - 1))
    elseif matrixtype == :banded
        matrix = BandedMatrix{Tv}(Zeros(n, n), (2 * nspec - 1, 2 * nspec - 1))
        # elseif matrixtype==:multidiagonal
        #     system.matrix=mdzeros(Tv,n,n,[-1,0,1]; blocksize=nspec)
    else # :sparse
        if num_partitions(system.grid) == 1
            matrix = ExtendableSparseMatrixCSC{Tv, Tm}(n, n)
        else
            matrix = MTExtendableSparseMatrixCSC{Tv, Tm}(n, n, num_partitions(system.grid))
        end
    end

    solution = unknowns(system)
    residual = unknowns(system)
    update = unknowns(system)
    dudp = [unknowns(system) for i = 1:(system.num_parameters)]
    SystemState(system, data, solution, matrix, dudp, residual, update, nothing, nothing, zero(UInt64))
end

SystemState(system::AbstractSystem{Tv,Tc,Ti,Tm}; kwargs...) where {Tv,Tc,Ti,Tm} =SystemState(Tv, system; kwargs...) 


# function _eval_and_assemble_inactive_species(system::AbstractSystem, U, Uold, F) end

# function _eval_and_assemble_inactive_species(system::DenseSystem, U, Uold, F)
#     if system.species_homogeneous
#         return
#     end
#     for inode = 1:size(system.node_dof, 2)
#         for ispec = 1:size(system.node_dof, 1)
#             if !isnodespecies(system, ispec, inode)
#                 F[ispec, inode] += U[ispec, inode] - Uold[ispec, inode]
#                 idof = dof(F, ispec, inode)
#                 rawupdateindex!(system.matrix, +, 1.0, idof, idof)
#             end
#         end
#     end
# end
# xX
# function _initialize_inactive_dof!(U::AbstractMatrix, system::AbstractSystem) end

# function _initialize_inactive_dof!(U::DenseSolutionArray, system::DenseSystem)
#     if system.species_homogeneous
#         return
#     end
#     for inode = 1:size(system.node_dof, 2)
#         for ispec = 1:size(system.node_dof, 1)
#             if !isnodespecies(system, ispec, inode)
#                 U[ispec, inode] = 0
#             end
#         end
#     end
# end

