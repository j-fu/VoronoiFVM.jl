################################################
const default_umfpack_pivot_tolerance=SuiteSparse.UMFPACK.umf_ctrl[3+1]

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
    tol_absolute::Float64 = 1.0e-10

    """
    Tolerance (relative to the size of the first update):
    terminate if ``\\Delta u_i/\\Delta u_1<`` `tol_relative`.
    """
    tol_relative::Float64 = 1.0e-10

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
    max_iterations::Int32 = 100

    """
    Maximum number of reuses of lu factorization.
    It this value is 0, linear systems are solved by a sparse direct solver, 
    and it's LU factorization is called in every Newton step.
    Otherwise, a BICGstab iterative method is used for linear system solution with a 
    LU factorization as preconditioner which is updated only every `max_lureuse` Newton step.
    """
    max_lureuse::Int32 = 0

    """
    Maximum number of consecutive iterations within roundoff error tolerance
    The default effectively disables this criterion.
    """
    max_round::Int32 = 1000

    """
    Tolerance of iterative linear solver.
    """
    tol_linear::Float64 = 1.0e-4

    """
    Factorization kind for linear sytems (see ExtendableSparse.jl). 
    Default: Standard Julia LU Factorization (UMFPACK).
    """
    factorization::AbstractFactorization = LUFactorization()

    """
    Iterative solver if factorization is incomplete.
    Currently supported: :bicgstab, :cg
    """
    iteration::Symbol=:bicgstab
    
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
    Pivot tolerance for umfpack.
    """
    umfpack_pivot_tolerance::Float64 =  default_umfpack_pivot_tolerance

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
    log=false
end


"""
````
timesteps!(control,Δt; grow=1.0)
````

Modify control data such that the time steps are fixed to a
geometric sequence such that Δt_new=Δt_old*grow
"""
function fixed_timesteps!(control,Δt; grow=1.0)
    control.Δt=Δt
    control.Δt_max=Δt
    control.Δt_min=Δt
    control.Δt_grow=grow
    control.Δu_opt=floatmax()
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
const NewtonControl=SolverControl


