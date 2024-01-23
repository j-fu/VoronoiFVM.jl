abstract type AbstractStrategy end

"""
    VoronoiFVM.LinearSolverStrategy

An linear solver strategy provides the possibility to construct [`SolverControl`](@ref) objects as follows:
```
    SolverControl(strategy,sys;kwargs...)
```,
e.g.
```
    SolverControl(GMRESIteration(UMFPackFactorization(), EquationBlock()),sys;kwargs...)
```

A linear solver strategy combines a Krylov method  with an (incomplete) factorization
which by default is calculated from the linearization of the initial value of the
Newton iteration.

Currently available strategies are:
- [`DirectSolver`](@ref)
- [`CGIteration`](@ref)
- [`BICGstabIteration`](@ref)
- [`GMRESIteration`](@ref)


Notable LU Factorizations are:
- [`UMFPACKFactorization`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/#SuiteSparse.jl)
- [`KLUFactorization`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/#SuiteSparse.jl)
- [`SparspakFactorization`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/#Sparspak.jl), [`SparspakLU`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.SparspakLU)
- [`MKLPardisoLU`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.MKLPardisoLU)

Notable incomplete factorizations are:
- [`ILUZeroPreconditioner`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.ILUZeroPreconditioner)
- [`AMGPreconditioner`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.AMGPreconditioner),
- [`ILUTPrecondidtioner`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.ILUTPreconditioner)

"""
abstract type LinearSolverStrategy<:AbstractStrategy end


"""
    VoronoiFVM.BlockStrategy

Abstract supertype for various block preconditioning strategies.
"""
abstract type BlockStrategy<:AbstractStrategy end

"""
    NoBlock()

No blocking.
"""
struct NoBlock <: BlockStrategy end

"""
    EquationBlock()

Equation-wise blocking. Can be combined with any preconditioner resp. factorization including direct solvers.
In the moment, this requires a system with `unknown_storage=:dense`.
"""
struct EquationBlock <: BlockStrategy end

"""
    PointBlock()

Point-wise blocking. Currently only together with ILUZeroFactorization.
This requires a system with `unknown_storage=:dense`.
"""
struct PointBlock <: BlockStrategy end


"""
    DirectSolver(;factorization=UMFPACKFactorization())

LU Factorization solver.
"""
Base.@kwdef struct DirectSolver <: LinearSolverStrategy
    factorization::FactorizationStrategy = UMFPACKFactorization()
    blocking::NoBlock = NoBlock() # prevent ambiguity in constructor definition
end

DirectSolver(factorization::FactorizationStrategy; kwargs...) =
    DirectSolver(; factorization, kwargs...)

function VoronoiFVM.SolverControl(strat::DirectSolver, sys; kwargs...)
    SolverControl(; method_linear = strat.factorization, kwargs...)
end




"""
    GMRESIteration(;factorization=ILUZeroFactorization(), memory=20, restart=true)
    GMRESIteration(factorization; memory=20, restart=true)
    
GMRES Iteration from Krylov.jl via LinearSolve.jl.
"""
Base.@kwdef struct GMRESIteration <: LinearSolverStrategy
    memory::Int = 20
    restart::Bool = true
    factorization::FactorizationStrategy = UMFPACKFactorization()
    blocking::BlockStrategy = NoBlock()
end


GMRESIteration(factorization::FactorizationStrategy, blocking = NoBlock(); kwargs...) =
    GMRESIteration(; factorization, blocking, kwargs...)

function VoronoiFVM.SolverControl(strat::GMRESIteration, sys; kwargs...)
    SolverControl(;
        method_linear = KrylovJL_GMRES(
            gmres_restart = strat.memory,
            restart = strat.restart,
        ),
        precon_linear = factorizationstrategy(strat.factorization, strat.blocking, sys),
        kwargs...,
    )
end


"""
    CGIteration(;factorization=UMFPACKFactorization())
    CGIteration(factorization)
    
CG Iteration from Krylov.jl via LinearSolve.jl.
"""
Base.@kwdef struct CGIteration <: LinearSolverStrategy
    factorization::FactorizationStrategy = UMFPACKFactorization()
    blocking::BlockStrategy = NoBlock()
end

CGIteration(factorization::FactorizationStrategy, blocking = NoBlock(); kwargs...) =
    CGIteration(; factorization, blocking, kwargs...)


function VoronoiFVM.SolverControl(strat::CGIteration, sys; kwargs...)
    SolverControl(;
        method_linear = KrylovJL_CG(),
        precon_linear = factorizationstrategy(strat.factorization, strat.blocking, sys),
        kwargs...,
    )
end


"""
    BICGstabIteration(;factorization=UMFPACKFactorization())
    BICGstabIteration(factorization)
    
BICGstab Iteration from Krylov.jl via LinearSolve.jl.
"""
Base.@kwdef struct BICGstabIteration <: LinearSolverStrategy
    factorization::FactorizationStrategy = UMFPACKFactorization()
    blocking::BlockStrategy = NoBlock()
end

BICGstabIteration(factorization::FactorizationStrategy, blocking = NoBlock(); kwargs...) =
    BICGstabIteration(; factorization, blocking, kwargs...)


function VoronoiFVM.SolverControl(strat::BICGstabIteration, sys; kwargs...)
    SolverControl(;
        method_linear = KrylovJL_BICGSTAB(),
        precon_linear = factorizationstrategy(strat.factorization, strat.blocking, sys),
        kwargs...,
    )
end



"""
    factorizationstrategy(preconditioner, blockstratrgy, system)

Create a factorizations strategy from preconditioner and block information
"""
factorizationstrategy(p::FactorizationStrategy, ::NoBlock, sys) = p

function factorizationstrategy(strat::FactorizationStrategy, ::EquationBlock, sys)
    !isdensesystem(sys) ?
    error("Equation block preconditioner currently needs dense system") : nothing
    BlockPreconditioner(;
        partitioning = partitioning(sys),
        factorization = factorizationstrategy(strat, NoBlock(), sys),
    )
end

function factorizationstrategy(::ILUZeroPreconditioner, ::PointBlock, sys)
    !isdensesystem(sys) ?
    error("Point block preconditioner needs dense system") : nothing
    PointBlockILUZeroPreconditioner(; blocksize = num_species(sys))
end

VoronoiFVM.SolverControl(::AbstractStrategy, sys; kwargs...) = SolverControl(;kwargs...)
VoronoiFVM.SolverControl(::Nothing, sys; kwargs...) = SolverControl(;kwargs...)


################################################################
# These are needed to enable iterative solvers to work with dual numbers
Base.Float64(x::ForwardDiff.Dual) = value(x)
function Random.rand(
    rng::AbstractRNG,
    ::Random.SamplerType{ForwardDiff.Dual{T,V,N}},
) where {T,V,N}
    ForwardDiff.Dual{T,V,N}(rand(rng, V))
end


"""
    Make preconditioner constructors from methods
"""

function (method::LinearSolve.AbstractFactorization)(A)
    pr = LinearProblem(A, zeros(eltype(A), size(A, 1)))
    init(pr, method)
end

function (method::LinearSolve.SciMLLinearSolveAlgorithm)(A)
    pr = LinearProblem(SparseMatrixCSC(A), zeros(eltype(A), size(A, 1)))
    init(pr, method)
end

function (f::ExtendableSparse.AbstractFactorization)(A)
    factorize!(f, A)
end

function LinearAlgebra.ldiv!(u, cache::LinearSolve.LinearCache, b)
    cache.b=b
    sol = solve!(cache)
    copyto!(u, sol.u)
end

canonical_matrix(A)=A
canonical_matrix(A::ExtendableSparseMatrix)=SparseMatrixCSC(A)

function _solve_linear!(u, system, nlhistory, control, method_linear, A, b)
    if isnothing(system.linear_cache)
        Pl = control.precon_linear(canonical_matrix(A))
        nlhistory.nlu += 1
        p = LinearProblem(canonical_matrix(A), b)
        system.linear_cache = init(
            p,
            method_linear;
            abstol = control.abstol_linear,
            reltol = control.reltol_linear,
            verbose = doprint(control, 'l'),
            Pl,
        )
    else
        system.linear_cache.A=canonical_matrix(A)
        system.linear_cache.b=b
        if control.keepcurrent_linear
            nlhistory.nlu += 1
            system.linear_cache.Pl=control.precon_linear(canonical_matrix(A))
        end
    end

    try
        sol = LinearSolve.solve!(system.linear_cache)
        u .= sol.u
        nliniter = sol.iters
        nlhistory.nlin = sol.iters
    catch err
        if (control.handle_exceptions)
            _print_error(err, stacktrace(catch_backtrace()))
            throw(LinearSolverError())
        else
            rethrow(err)
        end
    end
end

#################################
"""
    SolverStrategies

!!! compat "Only available in 1.5"
    Please replace this functionality by the new strategy API in 1.6 as follows:
    ```                                                                                
    direct_umfpack() = DirectSolver(UMFPACKFactorization())                            
    gmres_umfpack() = GMRESIteration(UMFPACKFactorization())                           
    gmres_eqnblock_umfpack() = GMRESIteration(UMFPACKFactorization(), EquationBlock()) 
    gmres_iluzero() = GMRESIteration(ILUZeroPreconditioner())                          
    gmres_eqnblock_iluzero() = GMRESIteration(ILUZeroPreconditioner(), EquationBlock())
    gmres_pointblock_iluzero() = GMRESIteration(ILUZeroPreconditioner(), PointBlock()) 
    ```                                                                                
"""
module SolverStrategies
using VoronoiFVM
using VoronoiFVM: isdensesystem, FactorizationStrategy
using DocStringExtensions
using LinearSolve
using ExtendableSparse



direct_umfpack() = DirectSolver(UMFPACKFactorization())
gmres_umfpack() = GMRESIteration(UMFPACKFactorization())
gmres_eqnblock_umfpack() = GMRESIteration(UMFPACKFactorization(), EquationBlock())
gmres_iluzero() = GMRESIteration(ILUZeroPreconditioner())
gmres_eqnblock_iluzero() = GMRESIteration(ILUZeroPreconditioner(), EquationBlock())
gmres_pointblock_iluzero() = GMRESIteration(ILUZeroPreconditioner(), PointBlock())

export direct_umfpack
export gmres_umfpack
export gmres_eqnblock_umfpack
export gmres_iluzero
export gmres_eqnblock_iluzero
export gmres_pointblock_iluzero


end
