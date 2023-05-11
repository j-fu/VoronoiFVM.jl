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

abstr

direct()
direct(umfpack())
gmres(eqnblock_iluzero())

Base.@kwdefstruct gmres <:IterationStrategy
    memory::Int = 20
    restart::Bool = true
    factorization<:FactoriationStrategy = umfpack()
end

"""


module SolverStrategies
using VoronoiFVM
using VoronoiFVM: isdensesystem
using DocStringExtensions
using LinearSolve
using ExtendableSparse

############ Direct solvers

for (strat, fact) in (
    (:direct_umfpack, :UMFPACKFactorization),
    (:direct_sparspak, :SparspakFactorization),
    (:direct_klu, :KLUFactorization),
)
    @eval begin
        @doc $("```\n$strat()\n```\n Direct solver via $fact.\n") struct $strat <:
                                                                         VoronoiFVM.AbstractStrategy end
        VoronoiFVM.SolverControl(::$strat, sys; kwargs...) =
            SolverControl(; method_linear = $fact(), kwargs...)
        export $strat
    end
end


for (strat, fact) in (
    (:gmres_umfpack, :UMFPACKFactorization),
    (:gmres_sparspak, :SparspakFactorization),
    (:gmres_klu, :KLUFactorization),
    (:gmres_iluzero, :ILUZeroPreconditioner),
)
    @eval begin
        @doc $(
            """```\n$strat(;memory=20, restart=true)\n```\n GMRES iterative solver with $fact preconditioner calculated at initial Newton iterate.\n
            Parameters:\n
            - `memory` maximal size of Krylov space before restart\n
            - `restart` allow restart"""
        ) Base.@kwdef struct $strat <: VoronoiFVM.AbstractStrategy
            memory::Int = 20
            restart::Bool = true
        end
        VoronoiFVM.SolverControl(strat::$strat, sys; kwargs...) = SolverControl(;
            method_linear = KrylovJL_GMRES(
                gmres_restart = strat.memory,
                restart = strat.restart,
            ),
            precon_linear = $fact(),
            kwargs...,
        )
        export $strat
    end
end

for (strat, fact) in (
    (:gmres_eqnblock_umfpack, :UMFPACKFactorization),
    (:gmres_eqnblock_sparspak, :SparspakFactorization),
    (:gmres_eqnblock_klu, :KLUFactorization),
    (:gmres_eqnblock_iluzero, :ILUZeroPreconditioner),
)
    @eval begin
        @doc $(
            """```\n$strat(;memory=20, restart=true)\n```\n GMRES iterative solver with equation block $fact preconditioner calculated at initial Newton iterate.\n
            Parameters:\n
            - `memory` maximal size of Krylov space before restart\n
            - `restart` allow restart"""
        ) Base.@kwdef struct $strat <: VoronoiFVM.AbstractStrategy
            memory::Int = 20
            restart::Bool = true
        end
        function VoronoiFVM.SolverControl(strat::$strat, sys; kwargs...)
            !isdensesystem(sys) ?
            error("Equation block preconditioner currently needs dense system") : nothing
            SolverControl(;
                method_linear = KrylovJL_GMRES(
                    gmres_restart = strat.memory,
                    restart = strat.restart,
                ),
                precon_linear = BlockPreconditioner(;
                    partitioning = partitioning(sys),
                    factorization = $fact(),
                ),
                kwargs...,
            )
        end
        export $strat
    end
end

for (strat, iter, fact) in (
    (:cg_umfpack, :KrylovJL_CG, :UMFPACKFactorization),
    (:cg_sparspak, :KrylovJL_CG, :SparspakFactorization),
    (:cg_klu, :KrylovJL_CG, :KLUFactorization),
    (:cg_iluzero, :KrylovJL_CG, :ILUZeroPreconditioner),
    (:bicgstab_umfpack, :KrylovJL_BICGSTAB, :UMFPACKFactorization),
    (:bicgstab_sparspak, :KrylovJL_BICGSTAB, :SparspakFactorization),
    (:bicgstab_klu, :KrylovJL_BICGSTAB, :KLUFactorization),
    (:bicgstab_iluzero, :KrylovJL_BICGSTAB, :ILUZeroPreconditioner),
)
    @eval begin
        @doc $(
            """```\n$strat()\n```\n $iter iterative solver with $fact preconditioner calculated at initial Newton iterate.\n"""
        ) Base.@kwdef struct $strat <: VoronoiFVM.AbstractStrategy end
        VoronoiFVM.SolverControl(strat::$strat, sys; kwargs...) =
            SolverControl(; method_linear = :iter(), precon_linear = $fact(), kwargs...)
        export $strat
    end
end

for (strat, iter, fact) in (
    (:cg_eqnblock_umfpack, :KrylovJL_CG, :UMFPACKFactorization),
    (:cg_eqnblock_sparspak, :KrylovJL_CG, :SparspakFactorization),
    (:cg_eqnblock_klu, :KrylovJL_CG, :KLUFactorization),
    (:cg_eqnblock_iluzero, :KrylovJL_CG, :ILUZeroPreconditioner),
    (:bicgstab_eqnblock_umfpack, :KrylovJL_BICGSTAB, :UMFPACKFactorization),
    (:bicgstab_eqnblock_sparspak, :KrylovJL_BICGSTAB, :SparspakFactorization),
    (:bicgstab_eqnblock_klu, :KrylovJL_BICGSTAB, :KLUFactorization),
    (:bicgstab_eqnblock_iluzero, :KrylovJL_BICGSTAB, :ILUZeroPreconditioner),
)
    @eval begin
        @doc $(
            """```\n$strat()\n```\n $iter iterative solver with $fact preconditioner calculated at initial Newton iterate.\n"""
        ) Base.@kwdef struct $strat <: VoronoiFVM.AbstractStrategy end
        function VoronoiFVM.SolverControl(strat::$strat, sys; kwargs...)
            !isdensesystem(sys) ?
            error("Equation block preconditioner currently needs dense system") : nothing
            SolverControl(;
                method_linear = :iter(),
                precon_linear = BlockPreconditioner(;
                    partitioning = partitioning(sys),
                    factorization = $fact(),
                ),
                kwargs...,
            )
        end
        export $strat
    end
end





"""
    $(TYPEDEF)

Solve linear systems during Newton's method using the GMRES iterative
method from Krylov.jl using an point-block preconditioner 
based on an incomplete LU factorization the  of the Jacobian of the initial value
see as a pointwise block matrix, calculated with ILUZero.jl.

Only available for dense system.
$(TYPEDFIELDS)
"""
Base.@kwdef struct gmres_pointblock_iluzero <: VoronoiFVM.AbstractStrategy
    memory::Int = 20
    restart::Bool = true
end
function VoronoiFVM.SolverControl(strat::gmres_pointblock_iluzero, sys; kwargs...)
    !isdensesystem(sys) ? error("Point block preconditioner needs dense system") : nothing
    SolverControl(;
        method_linear = KrylovJL_GMRES(
            gmres_restart = strat.memory,
            restart = strat.restart,
        ),
        precon_linear = PointBlockILUZeroPreconditioner(; blocksize = num_species(sys)),
    )
end
export gmres_pointblock_iluzero



for (strat, iter) in (
    (:cg_pointblock_iluzero, :KrylovJL_CG),
    (:bicgstab_pointblock_iluzero, :KrylovJL_BICGSTAB),
)
    @eval begin
        @doc $(
            """```\n$strat()\n```\n $iter  iterative solver with point block ILUZero preconditioner calculated at initial Newton iterate.\n"""
        ) Base.@kwdef struct $strat <: VoronoiFVM.AbstractStrategy end
        function VoronoiFVM.SolverControl(strat::$strat, sys; kwargs...)
            !isdensesystem(sys) ? error("Point block preconditioner needs dense system") :
            nothing
            SolverControl(;
                method_linear = $iter(),
                precon_linear = PointBlockILUZeroPreconditioner(;
                    blocksize = num_species(sys),
                ),
            )
        end
        export $strat
    end
end



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


# export gmres_amg
# export gmres_eqnblock_amg

#= Idea:
define all strategies here and export them.
Add SolverControl methods in corresponding extensions.
Probably do this only for 1.9 as this is less work.
=#

end
