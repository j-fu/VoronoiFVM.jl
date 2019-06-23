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
    Tolerance of iterative linear solver
    """
    tol_linear::Float64

    """
    Verbosity flag
    """
    verbose::Bool      
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
    this.damp_initial=1.0
    this.damp_growth=1.2
    this.max_lureuse=0
    this.tol_linear=1.0e-4
    this.verbose=false
    this.max_iterations=100
    return this
end

