abstract type AbstractStrategy end

"""
    VoronoiFVM.LinearSolverStrategy

An linear solver strategy provides the possibility to construct [`SolverControl`](@ref) objects as follows:
```
    SolverControl(strategy,sys;kwargs...)
```,
e.g.
```
    SolverControl(GMRESIteration(UMFPackFactorization(), EquationBlock()),sys; kwargs...)
```

A linear solver strategy combines a Krylov method  with a preconditioner
which by default is calculated from the linearization of the initial value of the
Newton iteration. For coupled systems, a blocking strategy can be chosen. The [`EquationBlock`](@ref) strategy
calculates preconditioners or LU factorization  separately for each species equation and combines them
to a block Jacobi preconditioner.  The [`PointBlock`](@ref) strategy treats the linear system as consisting
of `nspecies x nspecies` blocks. 

Which is the best strategy, depends on the space dimension. The following is a rule of thumb for choosing strategies
- For 1D problems use direct solvers
- For 2D stationary problems, use direct solvers, for transient problems consider iterative solvers which 
  can take advantage of the diagonal dominance of the implicit timestep problem
- For 3D problems avoid direct solvers


Currently available strategies are:
- [`DirectSolver`](@ref)
- [`CGIteration`](@ref)
- [`BICGstabIteration`](@ref)
- [`GMRESIteration`](@ref)

Notable LU Factorizations/direct solvers are:
- [`UMFPACKFactorization`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/#SuiteSparse.jl)  (`using LinearSolve`)
- [`KLUFactorization`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/#SuiteSparse.jl) (`using LinearSolve`)
- [`SparspakFactorization`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/#Sparspak.jl)  (`using LinearSolve`), [`SparspakLU`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.SparspakLU) (`using ExtendableSparse,Sparspak`)
- [`MKLPardisoLU`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.MKLPardisoLU) (`using ExtendableSparse, Pardiso`), openmp parallel
- [`AMGSolver`](https://j-fu.github.io/AMGCLWrap.jl/stable/solvers/#AMGCLWrap.AMGSolver) (`using AMGCLWrap`), openmp parallel
- [`RLXSolver`](https://j-fu.github.io/AMGCLWrap.jl/stable/solvers/#AMGCLWrap.RLXSolver) (`using AMGCLWrap`), openmp parallel


Notable incomplete factorizations/preconditioners
- Incomplete LU factorizations written in Julia:
    - [`ILUZeroPreconditioner`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.ILUZeroPreconditioner)
    - [`ILUTPrecondidtioner`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.ILUTPreconditioner) (`using ExtendableSparse, IncompleteLU`)
- Algebraic multigrid written in Julia: (`using ExtendableSparse, AlgebraicMultigrid`)
    - [`RS_AMGPreconditioner`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.RS_AMGPreconditioner)
    - [`SA_AMGPreconditioner`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.SA_AMGPreconditioner)
- AMGCL based preconditioners (`using ExtendableSparse, AMGCLWrap`), openmp parallel
    - [`AMGCL_AMGPreconditioner`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.AMGCL_AMGPreconditioner)
    - [`AMGCL_RLXPreconditioner`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.AMGCL_RLXPreconditioner)

Blocking strategies are:
- [`NoBlock`](@ref)
- [`EquationBlock`](@ref)
- [`PointBlock`](@ref)

"""
abstract type LinearSolverStrategy <: AbstractStrategy end

"""
    VoronoiFVM.BlockStrategy

Abstract supertype for various block preconditioning strategies.
"""
abstract type BlockStrategy <: AbstractStrategy end

"""
    NoBlock()

No blocking.
"""
struct NoBlock <: BlockStrategy end

"""
    EquationBlock()

Equation-wise blocking. Can be combined with any preconditioner resp. factorization including direct solvers.
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

DirectSolver(factorization::FactorizationStrategy; kwargs...) = DirectSolver(; factorization, kwargs...)

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

function GMRESIteration(factorization::FactorizationStrategy, blocking = NoBlock(); kwargs...)
    GMRESIteration(; factorization, blocking, kwargs...)
end

function VoronoiFVM.SolverControl(strat::GMRESIteration, sys; kwargs...)
    SolverControl(;
                  method_linear = KrylovJL_GMRES(; gmres_restart = strat.memory,
                                                 restart = strat.restart),
                  precon_linear = factorizationstrategy(strat.factorization, strat.blocking, sys),
                  kwargs...,)
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

function CGIteration(factorization::FactorizationStrategy, blocking = NoBlock(); kwargs...)
    CGIteration(; factorization, blocking, kwargs...)
end

function VoronoiFVM.SolverControl(strat::CGIteration, sys; kwargs...)
    SolverControl(;
                  method_linear = KrylovJL_CG(),
                  precon_linear = factorizationstrategy(strat.factorization, strat.blocking, sys),
                  kwargs...,)
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

function BICGstabIteration(factorization::FactorizationStrategy, blocking = NoBlock(); kwargs...)
    BICGstabIteration(; factorization, blocking, kwargs...)
end

function VoronoiFVM.SolverControl(strat::BICGstabIteration, sys; kwargs...)
    SolverControl(;
                  method_linear = KrylovJL_BICGSTAB(),
                  precon_linear = factorizationstrategy(strat.factorization, strat.blocking, sys),
                  kwargs...,)
end

"""
    factorizationstrategy(preconditioner, blockstratrgy, system)

Create a factorizations strategy from preconditioner and block information
"""
factorizationstrategy(p::FactorizationStrategy, ::NoBlock, sys) = p

function factorizationstrategy(strat::FactorizationStrategy, ::EquationBlock, sys)
    BlockPreconditioner(;
                        partitioning = partitioning(sys, Equationwise()),
                        factorization = factorizationstrategy(strat, NoBlock(), sys),)
end

function factorizationstrategy(::ExtendableSparse.ILUZeroPreconditioner, ::PointBlock, sys)
    !isdensesystem(sys) ?
    error("Point block preconditioner needs dense system") : nothing
    PointBlockILUZeroPreconditioner(; blocksize = num_species(sys))
end

VoronoiFVM.SolverControl(::AbstractStrategy, sys; kwargs...) = SolverControl(; kwargs...)
VoronoiFVM.SolverControl(::Nothing, sys; kwargs...) = SolverControl(; kwargs...)

################################################################
# These are needed to enable iterative solvers to work with dual numbers
Base.Float64(x::ForwardDiff.Dual) = value(x)
function Random.rand(rng::AbstractRNG,
                     ::Random.SamplerType{ForwardDiff.Dual{T, V, N}}) where {T, V, N}
    ForwardDiff.Dual{T, V, N}(rand(rng, V))
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
    cache.b = b
    sol = solve!(cache)
    copyto!(u, sol.u)
end

canonical_matrix(A) = A
canonical_matrix(A::AbstractExtendableSparseMatrixCSC) = SparseMatrixCSC(A)

function _solve_linear!(u, system, nlhistory, control, method_linear, A, b)
    if isnothing(system.linear_cache)
        Pl = control.precon_linear(canonical_matrix(A))
        nlhistory.nlu += 1
        p = LinearProblem(canonical_matrix(A), b)
        system.linear_cache = init(p,
                                   method_linear;
                                   abstol = control.abstol_linear,
                                   reltol = control.reltol_linear,
                                   verbose = doprint(control, 'l'),
                                   Pl,)
    else
        system.linear_cache.A = canonical_matrix(A)
        system.linear_cache.b = b
        if control.keepcurrent_linear
            nlhistory.nlu += 1
            system.linear_cache.Pl = control.precon_linear(canonical_matrix(A))
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
