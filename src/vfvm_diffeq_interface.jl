#
# Interface to VoronoiFVMDiffEq.jl
#
# For VoronoiFVM v0.18, v0.19 we allow breaking changes on
# this part of the API in patch revisions.
#
"""
    $(TYPEDEF)

History information for DiffEqHistory
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct DiffEqHistory
    """ number of combined jacobi/rhs evaluations"""
    nd = 0
    """ number of combined jacobi evaluations"""
    njac = 0
    """ number of rhs evaluations"""
    nf = 0
end

details(h::DiffEqHistory) = (nd = h.nd, njac = h.njac, nf = h.nf)
Base.summary(h::DiffEqHistory) = details(h)

"""
    $(SIGNATURES)

Evaluate functiaon and Jacobian at u if they have not been evaluated before at u.
See https://github.com/SciML/DifferentialEquations.jl/issues/521 for discussion of another way to do this.
"""
function _eval_res_jac!(sys, u, t)
    uhash = hash(u)
    if uhash != sys.uhash
        ur = reshape(u, sys)
        eval_and_assemble(sys, ur, ur, sys.residual, value(t), Inf, 0.0, zeros(0))
        sys.uhash = uhash
        sys.history.nd += 1
    end
end

"""
$(SIGNATURES)

Interpret the  discrete problem as an ODE/DAE problem. Provide the 
rhs function for DifferentialEquations.jl.
"""
function eval_rhs!(du, u, sys, t)
    _eval_res_jac!(sys, u, t)
    du .= -vec(sys.residual)
    sys.history.nf += 1
    nothing
end

"""
$(SIGNATURES)

Interpret the  discrete problem as an ODE/DAE problem. Provide the 
jacobi matrix calculation function for DifferentialEquations.jl.
"""
function eval_jacobian!(J, u, sys, t)
    _eval_res_jac!(sys, u, t)
    # Need to implement broadcast for ExtendableSparse.
    J .= -sys.matrix.cscmatrix
    sys.history.njac += 1
    nothing
end

"""
$(SIGNATURES)

Calculate the mass matrix for use with DifferentialEquations.jl.
Return a Diagonal matrix if it occurs to be diagonal, otherwise return a SparseMatrixCSC.
"""
function mass_matrix(system::AbstractSystem{Tv, Tc, Ti, Tm}) where {Tv, Tc, Ti, Tm}
    physics = system.physics
    data = physics.data
    node = Node(system)
    bnode = BNode(system)
    nspecies = num_species(system)
    ndof = num_dof(system)

    stor_eval = ResJacEvaluator(physics, :storage, zeros(Tv, nspecies), node, nspecies)
    bstor_eval = ResJacEvaluator(physics, :bstorage, zeros(Tv, nspecies), node, nspecies)

    U = unknowns(system; inival = 0)
    M = ExtendableSparseMatrix{Tv, Tm}(ndof, ndof)


    asm_res(idof, ispec) = nothing
    asm_param(idof, ispec, iparam) = nothing

    for item in nodebatch(system.assembly_data)
        for inode in noderange(system.assembly_data,item)
            _fill!(node,system.assembly_data,inode,item)
            @views evaluate!(stor_eval, U[:, node.index])
            jac_stor = jac(stor_eval)
            asm_jac(idof, jdof, ispec, jspec) = _addnz(M, idof, jdof, jac_stor[ispec, jspec], node.fac)
            assemble_res_jac(node, system, asm_res, asm_jac, asm_param)
        end
    end
    
    if isnontrivial(bstor_eval)
        for item in nodebatch(system.boundary_assembly_data)
            for ibnode in noderange(system.boundary_assembly_data,item)
                _fill!(bnode,system.boundary_assembly_data,ibnode,item)
                K = bnode.index
                @views evaluate!(bstor_eval, U[:, K])
                jac_bstor = jac(bstor_eval)
                asm_jac(idof, jdof, ispec, jspec) = _addnz(M, idof, jdof, jac_bstor[ispec, jspec], bnode.fac)
                assemble_res_jac(node, system, asm_res, asm_jac, asm_param)
            end
        end
    end
    Mcsc = SparseMatrixCSC(M)
    isdiag(Mcsc) ? Diagonal([Mcsc[i, i] for i = 1:ndof]) : Mcsc
end

"""
$(SIGNATURES)
Prepare system for use with VoronoiFVMDiffEq.

- `jacval`: value at which to evaluate jacobian to obtatin prototype
- `tjac`: time moment for jacobian
 
Returns a prototype for the jacobian.
"""
function prepare_diffeq!(sys, jacval, tjac)
    sys.history = DiffEqHistory()
    _complete!(sys; create_newtonvectors = true)
    _eval_res_jac!(sys, jacval, tjac)
    flush!(sys.matrix)
    sys.matrix.cscmatrix
end
