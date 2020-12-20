#
# Interface to DifferentialEquations.jl
#

#
# Evaluate function and Jacobian at u if they have not
# been evaluated before at u.
#
# See https://github.com/SciML/DifferentialEquations.jl/issues/521
# for discussion of another way to do this
#
function _eval_res_jac!(sys,u)
    uhash=hash(u)
    if uhash!=sys.uhash
        ur=reshape(u,sys)
        eval_and_assemble(sys,ur,ur,sys.residual,Inf)
        sys.uhash=uhash
    end
end

"""
$(SIGNATURES)

Assume the discrete problem is an ODE problem. Provide the 
rhs function for DifferentialEquations.jl.
"""
function eval_rhs!(du, u, sys,t)
    _eval_res_jac!(sys,u)
    du.=-vec(sys.residual)
    nothing
end

"""
$(SIGNATURES)

Assume the discrete problem is an ODE problem. Provide the 
jacobi matrix calculation function for DifferentialEquations.jl.
"""
function eval_jacobian!(J, u, sys,t)
    _eval_res_jac!(sys,u)
    J.=-sys.matrix
    nothing
end

"""
$(SIGNATURES)

Calculate the mass matrix for use with DifferentialEquations.jl.
"""
function mass_matrix(this::AbstractSystem)
    grid=this.grid
    nodeparams=(node,)
    coord=grid[Coordinates]
    
    geom=grid[CellGeometries][1]
    cellnodes=grid[CellNodes]
    mmatrix=zeros(size(coord,2))

    # Eventually take into account storage term here.
    for icell=1:num_cells(grid)
        for inode=1:num_nodes(geom)
            mmatrix[cellnodes[inode,icell]]+=this.cellnodefactors[inode,icell]
        end
    end
    return Diagonal(mmatrix)
end

"""
$(SIGNATURES)

Complete the system and provide the jacobi matrix as prototype
for the Jacobian.
"""
function jac_prototype(sys::AbstractSystem)
    _complete!(sys, create_newtonvectors=true)
    sys.matrix
end

