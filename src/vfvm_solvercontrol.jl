################################################

"""
$(TYPEDEF)

Solver control parameters for time stepping, embedding, Newton method control.
All field names can be used as keyword arguments for [`solve(system::VoronoiFVM.AbstractSystem; kwargs...)`](@ref)

Newton's method solves ``F(u)=0`` by the iterative procedure ``u_{i+1}=u_{i} - d_i F'(u_i)^{-1}F(u_i)``
starting with some inital value ``u_0``, where ``d_i`` is a damping parameter.


$(TYPEDFIELDS)
"""
@with_kw mutable struct SolverControl
    """
    Tolerance (in terms of norm of Newton update):  
    terminate if ``\\Delta u_i=||u_{i+1}-u_i||_\\infty <`` `tol_absolute`.
    """
    abstol::Float64 = 1.0e-10

    """
    Tolerance (relative to the size of the first update):
    terminate if ``\\Delta u_i/\\Delta u_1<`` `tol_relative`.
    """
    reltol::Float64 = 1.0e-10

    """
    Tolerance for roundoff error detection:
    terminate if   ``|\\;||u_{i+1}||_1 - ||u_{i}||_1\\;|/ ||u_{i}||_1<`` `tol_round` occured `max_round` times in a row.
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
    Maximum number of iterations.
    """
    maxiters::Int = 100

    """
    Maximum number of reuses of lu factorization.
    It this value is 0, linear systems are solved by a sparse direct solver, 
    and it's LU factorization is called in every Newton step.
    Otherwise, a BICGstab iterative method is used for linear system solution with a 
    LU factorization as preconditioner which is updated only every `max_lureuse` Newton step.
    """
    max_lureuse::Int = 0

    """
    Maximum number of consecutive iterations within roundoff error tolerance
    The default effectively disables this criterion.
    """
    max_round::Int = 1000

    """
    Solver kind for linear systems (see LinearSolve.jl).
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
    Preconditioner for linear systems
    """
    precon_linear::Union{Nothing, Symbol, Function} = nothing

    """
    Verbosity flag.
    """
    verbose::Bool = false

    """
    Handle exceptions during transient solver and parameter embedding. 
    If `true`, exceptions in Newton solves are catched, the embedding resp. time step is lowered, 
    and solution is retried.  
    """
    handle_exceptions::Bool = false

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
    Edge parameter cutoff for rectangular triangles.
    """
    edge_cutoff::Float64 = 0.0

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

    tol_absolute::Union{Float64, Nothing} = nothing
    tol_relative::Union{Float64, Nothing} = nothing
    damp::Union{Float64, Nothing} = nothing
    damp_grow::Union{Float64, Nothing} = nothing
    max_iterations::Union{Float64, Nothing} = nothing
    tol_linear::Union{Float64, Nothing} = nothing
end

const key_replacements = Dict(:tol_absolute => :abstol,
                              :tol_relative => :reltol,
                              :damp => :damp_initial,
                              :damp_grow => :damp_growth,
                              :max_iterations => :maxiters,
                              :tol_linear => :reltol_linear)

function fix_deprecations!(control)
    # compatibility to names in SolverControl which cannot be deprecated.
    for key ∈ keys(key_replacements)
        value = getproperty(control, key)
        if !isnothing(value)
            @warn "deprecated SolverControl entry $(key). Please replace by $(key_replacements[key])."
            setproperty!(control, key_replacements[key], value)
            setproperty!(control, key, nothing)
        end
    end
end

"""
````
timesteps!(control,Δt; grow=1.0)
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
