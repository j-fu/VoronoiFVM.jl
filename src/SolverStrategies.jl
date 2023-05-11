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
