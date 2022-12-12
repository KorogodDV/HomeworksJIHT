module CTSolverStability

using StructArrays
using ForwardDiff
using LinearAlgebra

export VanDerWaalsMixture
export VanDerWaalsComponent

include("types.jl")
include("nlsolve.jl")
include("contants.jl")
include("initials.jl")
include("solve_cubic.jl")
include("stability.jl")
include("vanderwaals.jl")

greet() = print("Hello World!")

end # module CTSolverStability
