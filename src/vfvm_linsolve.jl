abstract type AbstractStrategy end

"""
    VoronoiFVM.IterationStrategy

An iteration strategy provides the possibility to construct SolverControl objects as follows:
```
    SolverControl(strategy,sys;kwargs...)
```,
e.g.
```
    SolverControl(GMRESIteration(UMFPackFactorization()),sys;kwargs...)
```

An iteration strategy combines a Krylov method  with an (incomplete) factorization
which by default is calculated from the linearization of the initial value of the
Newton iteration.

Notable LU Factorizations are:
[`UMFPACKFactorization`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/#SuiteSparse.jl),
[`KLUFactorization`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/#SuiteSparse.jl),
[`SparspakFactorization`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/#Sparspak.jl),
[`SparspakLU`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.SparspakLU),
[`MKLPardisoLU`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.MKLPardisoLU)

Notable incomplete factorizations are:
[`ILUZeroPreconditioner`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.ILUZeroPreconditioner)
[`AMGPreconditioner`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.AMGPreconditioner),
[`ILUTPrecondidtioner`](https://j-fu.github.io/ExtendableSparse.jl/stable/iter/#ExtendableSparse.ILUTPreconditioner)

"""
abstract type IterationStrategy<:AbstractStrategy end


"""
    VoronoiFVM.BlockStrategy

Abstract supertype for various block preconditioning strategies.
"""
abstract type BlockStrategy<:AbstractStrategy end


struct NoBlock <: BlockStrategy end
struct EquationBlock <: BlockStrategy end
struct PointBlock <: BlockStrategy end



Base.@kwdef struct DirectSolver <: IterationStrategy
    factorization::FactorizationStrategy = UMFPACKFactorization()
    x::Nothing = nothing # prevent ambiguity in constructor definition
end

DirectSolver(factorization::FactorizationStrategy; kwargs...) =
    DirectSolver(; factorization, kwargs...)

function VoronoiFVM.SolverControl(strat::DirectSolver, sys; kwargs...)
    SolverControl(; method_linear = strat.factorization, kwargs...)
end




Base.@kwdef struct GMRESIteration <: IterationStrategy
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


Base.@kwdef struct CGIteration <: IterationStrategy
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


mutable struct FactorizationPreconditioner{C}
    cache::C
end

"""
   factorize(A,method::LinearSolve.AbstractFactorization)

Calculate an LU factorization of A using one of the methods available in LinearSolve.jl. 
"""
function LinearAlgebra.factorize(A, method::LinearSolve.AbstractFactorization)
    pr = LinearProblem(A, zeros(eltype(A), size(A, 1)))
    p = FactorizationPreconditioner(init(pr, method))
    p
end

function LinearAlgebra.factorize(A, method::LinearSolve.SciMLLinearSolveAlgorithm)
    pr = LinearProblem(SparseMatrixCSC(A), zeros(eltype(A), size(A, 1)))
    p = FactorizationPreconditioner(init(pr, method))
    p
end


function (f::LinearSolve.AbstractFactorization)(A)
    factorize(A, f)
end

function (f::LinearSolve.SciMLLinearSolveAlgorithm)(A)
    factorize(SparseMatrixCSC(A), f)
end

function (f::ExtendableSparse.AbstractFactorization)(A)
    factorize!(f, A)
end


function LinearAlgebra.ldiv!(u, p::FactorizationPreconditioner, b)
    p.cache = LinearSolve.set_b(p.cache, b)
    sol = solve(p.cache)
    p.cache = sol.cache
    copyto!(u, sol.u)
end

function _solve_linear!(u, system, nlhistory, control, method_linear, A, b)
    if isnothing(system.linear_cache)
        Pl = control.precon_linear(SparseMatrixCSC(A))
        nlhistory.nlu += 1
        p = LinearProblem(SparseMatrixCSC(A), b)
        system.linear_cache = init(
            p,
            method_linear;
            abstol = control.abstol_linear,
            reltol = control.reltol_linear,
            verbose = doprint(control, 'l'),
            Pl,
        )
    else
        system.linear_cache = LinearSolve.set_A(system.linear_cache, SparseMatrixCSC(A))
        system.linear_cache = LinearSolve.set_b(system.linear_cache, b)
        if control.keepcurrent_linear
            Pl = control.precon_linear(SparseMatrixCSC(A))
            nlhistory.nlu += 1
            system.linear_cache =
                LinearSolve.set_prec(system.linear_cache, Pl, LinearSolve.Identity())
        end
    end

    try
        sol = LinearSolve.solve(system.linear_cache)
        system.linear_cache = sol.cache
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
