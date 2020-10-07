################################################
const default_umfpack_pivot_tolerance=SuiteSparse.UMFPACK.umf_ctrl[3+1]

"""
$(TYPEDEF)

Control parameters for Newton method.

Newton's method solves ``F(u)=0`` by the iterative procedure ``u_{i+1}=u_{i} - d_i F'(u_i)^{-1}F(u_i)``
starting with some inital value ``u_0``, where ``d_i`` is the damping. 

$(TYPEDFIELDS)
"""
mutable struct NewtonControl

    """
    Tolerance (in terms of norm of Newton update):  
    terminate if ``\\Delta_i=||u_{i+1}-u_i||_\\infty <`` `tol_absolute`.

    Default value: `1.0e-10`.
    """
    tol_absolute::Float64

    """
    Tolerance (relative to the first residual):
    terminate if ``\\Delta_i/\\Delta_0<`` `tol_relative`.

    Default value: `1.0e-10`.
    """
    tol_relative::Float64

    """
    Tolerance for roundoff error detection:
    terminate if   ``|\\;||u_{i+1}||_1 - ||u_{i}||_1\\;|/ ||u_{i}||_1<`` `tol_round` occured `max_round` times in a row.

    Default value: `1.0e-10`.
    """
    tol_round::Float64

    """
    Tolerance for monotonicity test:
    terminat with error if ``\\Delta_i/\\Delta_{i-1}>`` `1/tol_mono`.

    Default value: 1.0e-3
    """
    tol_mono::Float64

    """
    Initial damping parameter ``d_0``.

    Default value: `1.0`.

    To handle convergence problems, set this to a value less than 1.
    """
    damp_initial::Float64

    """
    Damping parameter growth factor: ``d_{i+1}=\\min(d_i\\cdot`` `max_growth` ``,1)``

    Default value: `1.2`

    Generally it should be set to a value between 1 and 2.
    """
    damp_growth::Float64

    """
    Maximum number of iterations.

    Default value: 100
    """
    max_iterations::Int32


    """
    Maximum number of reuses of lu factorization.
    It this value is 0, linear systems are solved by a sparse direct solver, and it's LU factorization is called in every Newton step.

    Otherwise, a BICGstab iterative method is used for linear system solution with an LU factorization as preconditioner which is updated only every `max_lureuse` Newton step.

    Default value: 0.
    """
    max_lureuse::Int32

    """
    Maximum number of consecutive iterations within roundoff error tolerance

    Default value: 1000 (effectively switching of this criterion).
    """
    max_round::Int32

    """
    Tolerance of iterative linear solver.

    Default value: `1.0e-4`.
    """
    tol_linear::Float64

    """
    Verbosity flag

    Default value: `false`.
    """
    verbose::Bool      

    """
    Handle exceptions in [`embed!`](@ref) and  [`evolve!`](@ref) methods. If `true`, exceptions in Newton solves are catched, the embedding resp. time step is lowered, and solution is retried.

    Default value: `false`
    """
    handle_exceptions::Bool
    
    """
    Initial parameter step for [`embed!`](@ref) method.

    Default value: `1.0`
    """
    Δp::Float64

    """
    Maximal parameter step for [`embed!`](@ref) method

    Default value: `1.0`
    """
    Δp_max::Float64
    
    """
    Minimal parameter step for [`embed!`](@ref) method.

    Default value: `1.0e-3`
    """
    Δp_min::Float64

    """
    Time step for [`evolve!`](@ref) method.
    
    Default value: 0.1
    """
    Δt::Float64

    """
    Maximal time step for [`evolve!`](@ref) method.
    
    Default value: 1
    """
    Δt_max::Float64

    """
    Minimal time step for [`evolve!`](@ref) method.
    
    Default value: `1.0e-3`
    """
    Δt_min::Float64

    """
    Maximal step size growth for  [`evolve!`](@ref) method.
    
    Default: 1.2
    """
    Δt_grow::Float64

    """
    Optimal size of update for  [`evolve!`](@ref) method.
    The algorithm keeps this value approximately constant.
    
    Default: 0.1
    """
    Δu_opt::Float64

    """
    Edge cutoff for rectangular triangles.

    Default value: `0.0`.
    """
    edge_cutoff::Float64

    """
    Pivot tolerance for umfpack

    Default value: $(default_umfpack_pivot_tolerance).
    """
    umfpack_pivot_tolerance::Float64

    
    NewtonControl()=NewtonControl(new())
end

################################################
"""
$(TYPEDSIGNATURES)
    
Default constructor
"""
function NewtonControl(this)
    this.tol_absolute=1.0e-10
    this.tol_relative=1.0e-10
    this.tol_round=1.0e-10
    this.tol_mono=1.0e-3
    this.damp_initial=1.0
    this.damp_growth=1.2
    this.max_lureuse=0
    this.tol_linear=1.0e-4
    this.max_round=1000
    this.verbose=false
    this.max_iterations=100
    this.Δp=1
    this.Δp_max=1
    this.Δp_min=1.0e-3
    this.Δt=0.1
    this.Δt_max=1
    this.Δt_min=1.0e-3
    this.Δt_grow=1.2
    this.Δu_opt=0.1
    this.handle_exceptions=false
    this.edge_cutoff=0.0
    this.umfpack_pivot_tolerance=default_umfpack_pivot_tolerance
    return this
end

function Base.show(io::IO, this::NewtonControl)
    for name in fieldnames(typeof(this))
        @printf("%16s = ",name)
        println(io,getfield(this,name))
    end
end
