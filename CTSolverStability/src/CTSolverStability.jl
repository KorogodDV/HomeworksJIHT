module CTSolverStability

using ForwardDiff

include("types.jl")
include("nlsolve.jl")
include("stability.jl")

greet() = print("Hello World!")

end # module CTSolverStability
