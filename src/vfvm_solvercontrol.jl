################################################

"""
    SolverControl
    SolverControl(;kwargs...)
    SolverControl(strategy, sys; kwargs...)

Solver control parameter for time stepping, embedding, Newton method and linear solver control.
All field names can be used as keyword arguments for [`solve(system::VoronoiFVM.AbstractSystem; kwargs...)`](@ref)

Newton's method solves ``F(u)=0`` by the iterative procedure ``u_{i+1}=u_{i} - d_i F'(u_i)^{-1}F(u_i)``
starting with some initial value ``u_0``, where ``d_i`` is a damping parameter.

Preset linear solver strategies are available from the submodule [`VoronoiFVM.SolverStrategies`](@ref).

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SolverControl
    """
    Verbosity control. A collection of output categories is given in a string composed of the
    following letters:
    -  a: allocation warnings
    -  d: deprecation warnings
    -  e: time/parameter evolution log
    -  n: newton solver log
    -  l: linear solver log
    Alternatively, a Bool value can be given, resulting in
    - true: "neda"
    - false: "da"
    Switch off all output including deprecation warnings via `verbose=""`.
    In the output, corresponding messages are marked e.g. via '[n]', `[a]` etc. (besides of '[l]')
    """
    verbose::Union{Bool, String} = false

    """
    Tolerance (in terms of norm of Newton update):  
    terminate if ``\\Delta u_i=||u_{i+1}-u_i||_\\infty <`` `abstol`.
    """
    abstol::Float64 = 1.0e-10

    """
    Tolerance (relative to the size of the first update):
    terminate if ``\\Delta u_i/\\Delta u_1<`` `reltol`.
    """
    reltol::Float64 = 1.0e-10

    """
    Maximum number of newton iterations.
    """
    maxiters::Int = 100

    """
    Tolerance for roundoff error detection:
    terminate if   ``|\\;||u_{i+1}||_1 - ||u_{i}||_1\\;|/ ||u_{i}||_1<`` `tol_round` occurred `max_round` times in a row.
    """
    tol_round::Float64 = 1.0e-10

    """
    Tolerance for monotonicity test:
    terminate with error if ``\\Delta u_i/\\Delta u_{i-1}>`` `1/tol_mono`.
    """
    tol_mono::Float64 = 1.0e-3

    """
    Initial damping parameter ``d_0``.
    To handle convergence problems, set this to a value less than 1.
    """
    damp_initial::Float64 = 1.0

    """
    Damping parameter growth factor: ``d_{i+1}=\\min(d_i\\cdot`` `max_growth` ``,1)``.
    It should be larger than 1.
    """
    damp_growth::Float64 = 1.2

    """
    Maximum number of consecutive iterations within roundoff error tolerance
    The default effectively disables this criterion.
    """
    max_round::Int = 1000

    """
    Calculation of Newton update norm
    """
    unorm::Function = (u) -> LinearAlgebra.norm(values(u), Inf) # norm for update calculation

    """
    Functional for roundoff error calculation
    """
    rnorm::Function = (u) -> LinearAlgebra.norm(values(u), 1)

    """
    Solver method for linear systems (see LinearSolve.jl). If given `nothing`, as default
    are chosen (for `Float64` calculations):
    - 1D:  `KLUFactorization()`
    - 2D:  `SparspakFactorization()`
    - 3D:  `UMFPACKFactorization()`
    `SparspakFactorization()` is the default choice for general number types.
    Users should experiment with what works best for their problem.
    """
    method_linear::Union{Nothing, LinearSolve.SciMLLinearSolveAlgorithm} = nothing

    """
        Relative tolerance of iterative linear solver.
    """
    reltol_linear::Float64 = 1.0e-4

    """
    Absolute tolerance of iterative linear solver.
    """
    abstol_linear::Float64 = 1.0e-8

    """
    Maximum number of iterations of linear solver
    """
    maxiters_linear::Int = 100

    """
    Constructor for preconditioner for linear systems.
    This should work as a function `precon_linear(A)` which
    takes an AbstractSparseMatrixCSC (e.g. an ExtendableSparseMatrix)
    and returns a preconditioner object in the sense of `LinearSolve.jl`, i.e. which
    has an `ldiv!(u,A,v)` method. Useful examples:
    - `ExtendableSparse.ILUZero`
    - `ExtendableSparse.Jacobi`
    """
    precon_linear::Union{Type, Function, ExtendableSparse.AbstractFactorization, LinearSolve.AbstractFactorization} = A -> Identity()

    """
    Update preconditioner in each Newton step ?
    """
    keepcurrent_linear::Bool = false

    """
    Initial parameter step for embedding.
    """
    Δp::Float64 = 1.0

    """
    Maximal parameter step size.
    """
    Δp_max::Float64 = 1.0

    """
    Minimal parameter step size.
    """
    Δp_min::Float64 = 1.0e-3

    """
    Maximal parameter step size growth. 
    """
    Δp_grow::Float64 = 1.0

    """
    Initial time step  size.
    """
    Δt::Float64 = 0.1

    """
    Maximal time step size. 
    """
    Δt_max::Float64 = 1.0

    """
    Minimal time step size.
    """
    Δt_min::Float64 = 1.0e-3

    """
    Maximal time step size growth.
    """
    Δt_grow::Float64 = 1.2

    """
    Optimal size of update for time stepping and embeding.
    The algorithm tries to keep the difference in norm between "old" and "new" 
    solutions  approximately at this value.
    """
    Δu_opt::Float64 = 0.1

    """
    Force first timestep.
    """
    force_first_step::Bool = false

    """
    Handle exceptions during transient solver and parameter embedding. 
    If `true`, exceptions in Newton solves are caught, the embedding resp. time step is lowered, 
    and solution is retried.  

    """
    handle_exceptions::Bool = false

    """
    Store all steps of transient/embedding problem:
    """
    store_all::Bool = true

    """
    Store transient/embedding solution in memory
    """
    in_memory::Bool = true

    """
       Record history
    """
    log = false

    """
    Edge parameter cutoff for rectangular triangles.
    """
    edge_cutoff::Float64 = 0.0

    """
    Function `pre(sol,t)` called before time/embedding step
    """
    pre::Function = function (sol, t) end

    """
    Function `post(sol,oldsol,t,Δt)` called after successful time/embedding step
    """
    post::Function = function (sol, oldsol, t, Δt) end

    """
    Function `sample(sol,t)` to be called for each `t in times[2:end]`
    """
    sample::Function = function (sol, t) end

    """
    Time step error estimator
    """
    delta::Function = (system, u, v, t, Δt) -> norm(system, u - v, Inf)

    # deprecated entries
    tol_absolute::Union{Float64, Nothing} = nothing
    tol_relative::Union{Float64, Nothing} = nothing
    damp::Union{Float64, Nothing} = nothing
    damp_grow::Union{Float64, Nothing} = nothing
    max_iterations::Union{Int, Nothing} = nothing
    tol_linear::Union{Float64, Nothing} = nothing
    max_lureuse::Union{Int, Nothing} = nothing
    mynorm::Union{Function, Nothing} = nothing
    myrnorm::Union{Function, Nothing} = nothing
end

doprint(s::String, a::Char) = contains(s, a)
const true_verbosity = "neda"
const false_verbosity = "da"
doprint(b::Bool, a::Char) = b ? doprint(true_verbosity, a) : doprint(false_verbosity, a)
doprint(c::SolverControl, a::Char) = doprint(c.verbose, a)

const key_replacements = Dict(:tol_absolute => :abstol,
                              :tol_relative => :reltol,
                              :damp => :damp_initial,
                              :damp_grow => :damp_growth,
                              :max_iterations => :maxiters,
                              :tol_linear => :reltol_linear,
                              :mynorm => :unorm,
                              :myrnorm => :rnorm,
                              :max_lureuse => nothing)

function fix_deprecations!(control)
    # compatibility to names in SolverControl which cannot be deprecated.
    for key ∈ keys(key_replacements)
        value = getproperty(control, key)
        if !isnothing(value)
            if !isnothing(key_replacements[key])
                if doprint(control, 'd')
                    @warn "[d]eprecated SolverControl entry '$(key)'. Please replace by '$(key_replacements[key])'."
                end
                setproperty!(control, key_replacements[key], value)
            else
                if doprint(control, 'd')
                    @warn "[d]eprecated SolverControl entry '$(key)' will be ignored."
                end
            end
            setproperty!(control, key, nothing)
        end
    end
end

"""
````
fixed_timesteps!(control,Δt; grow=1.0)
````

Modify control data such that the time steps are fixed to a
geometric sequence such that Δt_new=Δt_old*grow
"""
function fixed_timesteps!(control, Δt; grow = 1.0)
    control.Δt = Δt
    control.Δt_max = Δt
    control.Δt_min = Δt
    control.Δt_grow = grow
    control.Δu_opt = floatmax()
    control
end

# function Base.show(io::IO, this::SolverControl)
#     for name in fieldnames(typeof(this))
#         println(io,"$(lpad(name,20)) = $(getfield(this,name))")
#     end
# end

"""
    NewtonControl

Legacy name of SolverControl
"""
const NewtonControl = SolverControl



abstract type AbstractStrategy end

VoronoiFVM.SolverControl(::AbstractStrategy, sys; kwargs...) = SolverControl(;kwargs...)
VoronoiFVM.SolverControl(::Nothing, sys; kwargs...) = SolverControl(;kwargs...)


"""
    VoronoiFVM.SolverStrategies

This module contains a number of strategies which help to instantiate [`SolverControl`](@ref) objects
with solution methods and preconditioners/factorizations. General usage is
```
    SolverControl(strategy,sys;kwargs...)
```
See below for currently implemented strategies.
"""
module SolverStrategies
using VoronoiFVM
using DocStringExtensions
using LinearSolve
using ExtendableSparse

"""
    $(TYPEDEF)

Solve linear systems during Newton's method using the UMFPACK default sparse
solver of Julia.
"""
struct direct_umfpack <: VoronoiFVM.AbstractStrategy end
VoronoiFVM.SolverControl(::direct_umfpack, sys; kwargs...) =
    SolverControl(; method_linear = UMFPACKFactorization(), kwargs...)

"""
    $(TYPEDEF)

Solve linear systems during Newton's method using the GMRES iterative
method from Krylov.jl using an LU factorization of the Jacobian of the intial value
calculated with UMFPACK  as preconditioner.
"""
struct gmres_umfpack <: VoronoiFVM.AbstractStrategy end
VoronoiFVM.SolverControl(::gmres_umfpack, sys; kwargs...) = SolverControl(;
    method_linear = KrylovJL_GMRES(),
    precon_linear = UMFPACKFactorization(),
    kwargs...,
)

"""
    $(TYPEDEF)

Solve linear systems during Newton's method using the GMRES iterative
method from Krylov.jl using an equation-block preconditioner 
based on an LU factorization the equation blocks of the Jacobian of the initial value
calculated with UMFPACK  as preconditioner.
"""
struct gmres_eqnblock_umfpack <: VoronoiFVM.AbstractStrategy end
VoronoiFVM.SolverControl(::gmres_eqnblock_umfpack, sys; kwargs...) = SolverControl(;
    method_linear = KrylovJL_GMRES(),
    precon_linear = BlockPreconditioner(;
        partitioning = partitioning(sys),
        factorization = UMFPACKFactorization(),
    ),
        kwargs...,
)


"""
    $(TYPEDEF)

Solve linear systems during Newton's method using the GMRES iterative
method from Krylov.jl using an zero-fillin inclomplete LU factorization of the Jacobian of the initial
value calculated with ILUZero.jl.
"""
struct gmres_iluzero <: VoronoiFVM.AbstractStrategy end
VoronoiFVM.SolverControl(::gmres_iluzero, sys; kwargs...) = SolverControl(;
    method_linear = KrylovJL_GMRES(),
    precon_linear = ILUZeroPreconditioner(),
    kwargs...,
)

"""
    $(TYPEDEF)

Solve linear systems during Newton's method using the GMRES iterative
method from Krylov.jl using an equation-block preconditioner 
based on an incomplete LU factorization the equation blocks of the Jacobian of the initial value
calculated with ILUZero.jl.
"""
struct gmres_eqnblock_iluzero <: VoronoiFVM.AbstractStrategy end
VoronoiFVM.SolverControl(::gmres_eqnblock_iluzero, sys; kwargs...) = SolverControl(;
    method_linear = KrylovJL_GMRES(),
    precon_linear = BlockPreconditioner(;
        partitioning = partitioning(sys),
        factorization = ILUZeroPreconditioner(),
                   ),
        kwargs...
           )


"""
    $(TYPEDEF)

Solve linear systems during Newton's method using the GMRES iterative
method from Krylov.jl using an point-block preconditioner 
based on an incomplete LU factorization the  of the Jacobian of the initial value
see as a pointwise block matrix, calculated with ILUZero.jl.
"""
struct gmres_pointblock_iluzero <: VoronoiFVM.AbstractStrategy end
VoronoiFVM.SolverControl(::gmres_pointblock_iluzero, sys; kwargs...) = SolverControl(;
    method_linear = KrylovJL_GMRES(),
    precon_linear = PointBlockILUZeroPreconditioner(;
        blocksize = num_species(sys),
    ),
        kwargs...,
)

# struct gmres_amg <: VoronoiFVM.AbstractStrategy end
# SolverControl(::gmres_amg, sys; kwargs...) = SolverControl(;
#     method_linear = KrylovJL_GMRES(),
#     precon_linear = AMGPreconditioner(),
#     kwargs...,
# )

# struct gmres_eqnblock_amg <: VoronoiFVM.AbstractStrategy end
# SolverControl(::gmres_eqnblock_amg, sys; kwargs...) = SolverControl(;
#     method_linear = KrylovJL_GMRES(),
#     precon_linear = BlockPreconditioner(
#         partitioning = partitioning(sys),
#         factorization = AMGPreconditioner(),
#         kwargs...,
#     ),
# )


export direct_umfpack
export gmres_umfpack
export gmres_eqnblock_umfpack
export gmres_iluzero
export gmres_eqnblock_iluzero
export gmres_pointblock_iluzero

# export gmres_amg
# export gmres_eqnblock_amg

#= Idea:
define all strategies here and export them.
Add SolverControl methods in corresponding extensions.
Probably do this only for 1.9 as this is less work.
=#

end
