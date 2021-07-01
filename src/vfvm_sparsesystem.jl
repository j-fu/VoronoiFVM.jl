const SparseSystem = System{Tv,Ti,Tm,SparseMatrixCSC{Ti,Ti},SparseSolutionArray{Tv,Ti} } where {Tv, Ti, Tm}

"""
$(SIGNATURES)

Constructor for SparseSystem.
"""
SparseSystem(grid,physics::Physics; matrixindextype=Int64)=System(grid,physics,matrixindextype=matrixindextype,unknown_storage=:sparse)


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
    return SparseSolutionArray{Tu,Ti}(SparseMatrixCSC(system.node_dof.m,
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

