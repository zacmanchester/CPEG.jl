using Test
using CPEG
using BenchmarkTools

import CPEG

@testset "Atmosphere Tests" begin
    include("atmotest.jl")
end
