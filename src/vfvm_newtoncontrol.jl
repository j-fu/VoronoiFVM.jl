################################################
"""
$(TYPEDEF)

Control parameter structure for Newton method.

$(TYPEDFIELDS)
"""
mutable struct NewtonControl

    """
    Tolerance (in terms of norm of Newton update)
    """
    tol_absolute::Float64

    """
    Tolerance (relative to the first residual)
    """
    tol_relative::Float64

    """
    Tolerance for roundoff error detection
    """
    tol_round::Float64

    """
    Initial damping parameter
    """
    damp_initial::Float64

    """
    Damping parameter growth factor
    """
    damp_growth::Float64

    """
    Maximum number of iterations
    """
    max_iterations::Int32

    """
    Maximum number of reuses of lu factorization
    """
    max_lureuse::Int32

    """
    Maximum number of consecutive iterations within roundoff error tolerance
    """
    max_round::Int32

    
    """
    Tolerance of iterative linear solver
    """
    tol_linear::Float64

    """
    Verbosity flag
    """
    verbose::Bool      

    """
    Handle exceptions
    """
    handle_exceptions::Bool
    
    """
    Parameter step for embedding
    """
    Δp::Float64

    """
    Maximal parameter step
    """
    Δp_max::Float64
    
    """
    Minimal parameter step
    """
    Δp_min::Float64


    """
    Edge cutoff for rectangular triangles
    """
    edge_cutoff::Float64

    Δt::Float64
    Δt_max::Float64
    Δt_min::Float64
    Δt_grow::Float64
    Δu_opt::Float64
    
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
    
    return this
end

function Base.show(io::IO, this::NewtonControl)
    for name in fieldnames(typeof(this))
        @printf("%16s = ",name)
        println(io,getfield(this,name))
    end
end
