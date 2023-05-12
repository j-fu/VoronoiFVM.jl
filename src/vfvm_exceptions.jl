"""
$(TYPEDEF)

Exception thrown if Newton's method convergence fails.
"""
struct ConvergenceError <: Exception end

"""
$(TYPEDEF)

Exception thrown if error occurred during assembly (e.g. domain error)
"""
struct AssemblyError <: Exception end

"""
$(TYPEDEF)

Exception thrown if error occurred during factorization.
"""
struct LinearSolverError <: Exception end

"""
$(TYPEDEF)

Exception thrown if embedding fails
"""
struct EmbeddingError <: Exception
    msg::String
end

"""
Print error when catching exceptions
"""
function _print_error(err, st)
    println()
    println(err)
    nlines = 5
    for i = 1:min(nlines, length(st))
        line = @sprintf("%s", st[i])
        L = length(line)
        if L < 80
            println(line)
        else
            print(line[1:35])
            print(" ... ")
            println(line[(L-35):L])
        end
    end
    if length(st) > nlines
        println("...")
    end
    println()
end
