__precompile__(true)
module CPEG

# greet() = print("Hello World!")

using LinearAlgebra
using StaticArrays
using ForwardDiff
using SparseArrays
using SuiteSparse

include("atmosphere.jl")
include("scaling.jl")
include("gravity.jl")
include("aero_forces.jl")
include("vehicle.jl")
include("dynamics.jl")
include("post_process.jl")

end # module
