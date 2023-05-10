"""
    VoronoiFVM.SolverStrategies

This module contains a number of strategies which help to instantiate [`SolverControl`](@ref) objects
with pre-set solution methods and preconditioners/factorizations. General usage is: 
```
    SolverControl(strategy,sys;kwargs...)
```,
e.g.
```
    using VoronoiFVM.SolverStrategies
    SolverControl(gmres_umfpack(),sys;kwargs...)
```
See below for currently implemented strategies.

!!! note "Experimental"
    Please consider this feature experimental in version 1.5, i.e.
    1.6 may introduce breaking change for this functionality.
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
