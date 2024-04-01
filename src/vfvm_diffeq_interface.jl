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
rhs function for [`ODEFunction`](@ref).
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
jacobi matrix calculation function for [`ODEFunction`](@ref)
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

Calculate the mass matrix for use with [`ODEFunction`](@ref).
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
        for inode in noderange(system.assembly_data, item)
            _fill!(node, system.assembly_data, inode, item)
            @views evaluate!(stor_eval, U[:, node.index])
            jac_stor = jac(stor_eval)
            asm_jac(idof, jdof, ispec, jspec) = _addnz(M, idof, jdof, jac_stor[ispec, jspec], node.fac)
            assemble_res_jac(node, system, asm_res, asm_jac, asm_param)
        end
    end

    if isnontrivial(bstor_eval)
        for item in nodebatch(system.boundary_assembly_data)
            for ibnode in noderange(system.boundary_assembly_data, item)
                _fill!(bnode, system.boundary_assembly_data, ibnode, item)
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

###################################################################################################
# API

"""
     ODEFunction(system,inival=unknowns(system,inival=0),t0=0)
    
Create an [ODEPFunction](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems)
in [mass matrix form](https://diffeq.sciml.ai/stable/solvers/dae_solve/#OrdinaryDiffEq.jl-(Mass-Matrix))
to be handeled by ODE solvers from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

Parameters:
- `system`: A [`VoronoiFVM.System`](https://j-fu.github.io/VoronoiFVM.jl/stable/system/#VoronoiFVM.System-Tuple{ExtendableGrid})
- `jacval` (optional): Initial value. Default is a zero vector. Consider to  pass a stationary solution at time `tjac`.
- `tjac` (optional): tjac, Default: 0

The `jacval` and `tjac` are passed  for a first evaluation of the Jacobian, allowing to detect
the sparsity pattern which is passed to the solver.
"""
function SciMLBase.ODEFunction(sys::VoronoiFVM.AbstractSystem; jacval = unknowns(sys, 0), tjac = 0)
    SciMLBase.ODEFunction(eval_rhs!;
                          jac = eval_jacobian!,
                          jac_prototype = prepare_diffeq!(sys, jacval, tjac),
                          mass_matrix = mass_matrix(sys))
end

"""
    ODEProblem(system,inival,tspan,callback=SciMLBase.CallbackSet())
    
Create an [ODEProblem](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems)
in [mass matrix form](https://diffeq.sciml.ai/stable/solvers/dae_solve/#OrdinaryDiffEq.jl-(Mass-Matrix))
which can  be handeled by ODE solvers from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

Parameters:
- `system`: A [`VoronoiFVM.System`](https://j-fu.github.io/VoronoiFVM.jl/stable/system/#VoronoiFVM.System-Tuple{ExtendableGrid})
- `inival`: Initial value. Consider to  pass a stationary solution at `tspan[1]`.
- `tspan`: Time interval 
- `callback` : (optional) [callback](https://diffeq.sciml.ai/stable/features/callback_functions/#Using-Callbacks) for ODE solver 

The method returns an [ODEProblem](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems) which can be solved
by [solve()](https://diffeq.sciml.ai/stable/basics/common_solver_opts/).
"""
function SciMLBase.ODEProblem(sys::VoronoiFVM.AbstractSystem, inival, tspan, callback = SciMLBase.CallbackSet())
    odefunction = SciMLBase.ODEFunction(sys; jacval = inival, tjac = tspan[1])
    SciMLBase.ODEProblem(odefunction, vec(inival), tspan, sys, callback)
end

"""
    reshape(ode_solution, system, times=nothing)
Create a [`TransientSolution`](@ref) from the output of the ode solver which
reflects the species structure of the system ignored by the ODE solver.
Howvever the interpolation behind `reshaped_sol(t)` will be linear and ignores the possibility
of higher order interpolations with `ode_sol`.

If `times` is specified, the (possibly higher ordee) interpolated solution at the given moments of time will be returned.
"""
function Base.reshape(sol::AbstractDiffEqArray, sys::VoronoiFVM.AbstractSystem, times = nothing)
    if isnothing(times)
        TransientSolution([reshape(sol.u[i], sys) for i = 1:length(sol.u)], sol.t)
    else
        isol = sol(times)
        TransientSolution([reshape(isol[i], sys) for t = 1:length(times)], times)
    end
end
