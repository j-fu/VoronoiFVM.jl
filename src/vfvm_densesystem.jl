const DenseSystem = System{Tv,Ti,Tm,Matrix{Ti},Matrix{Tv}} where {Tv, Ti, Tm}


##################################################################
"""
$(SIGNATURES)

Constructor for DenseSystem.
"""
DenseSystem(grid,physics::Physics; matrixindextype=Int64)=System(grid,physics,matrixindextype=matrixindextype,unknown_storage=:dense)


##################################################################
"""
$(SIGNATURES)

Number of degrees of freedom for system.
"""
num_dof(system::DenseSystem)= length(system.node_dof)



##################################################################
"""
$(SIGNATURES)

Create a solution vector for dense system.
If inival is not specified, the entries of the returned vector are undefined.
"""
unknowns(system::DenseSystem{Tv};inival=undef) where Tv = unknowns(Tv,system,inival=inival)

##################################################################
"""
$(SIGNATURES)

Create a solution vector for dense system with elements of type `Tu`.
If inival is not specified, the entries of the returned vector are undefined.
"""
function unknowns(Tu::Type, sys::DenseSystem{Tv}; inival=undef) where Tv
    a=Array{Tu}(undef,size(sys.node_dof)...)
    if inival!=undef
        fill!(a,inival)
    end
    return a
end


Base.reshape(v::DenseSolutionArray{Tu},system::DenseSystem{Tv}) where {Tu,Tv}=v


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
function _eval_and_assemble_inactive_species(system::DenseSystem,U,Uold,F)
    if system.species_homogeneous
        return
    end
    for inode=1:size(system.node_dof,2)
        for ispec=1:size(system.node_dof,1)
            if !isdof(system,ispec,inode)
                F[ispec,inode]+= U[ispec,inode]-Uold[ispec,inode];
                idof=dof(F,ispec,inode)
                system.matrix[idof,idof]+=1.0
            end
        end
    end
end

#
# Initialize values in inactive dof for dense system
#
function _initialize_inactive_dof!(U::DenseSolutionArray{Tv},system::DenseSystem{Tv}) where {Tv}
    if system.species_homogeneous
        return
    end
    for inode=1:size(system.node_dof,2)
        for ispec=1:size(system.node_dof,1)
            if !isdof(system,ispec,inode)
                U[ispec,inode]=0
            end
        end
    end
end


"""
$(SIGNATURES)

Calculate norm, paying attention to species distribution over regions
"""
function LinearAlgebra.norm(system::DenseSystem,u,p)
    _initialize_inactive_dof!(u,system)
    norm(u,p)
end

