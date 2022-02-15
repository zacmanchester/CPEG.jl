using Test
using CPEG
using BenchmarkTools
using StaticArrays
using LinearAlgebra

import CPEG

@testset "Atmosphere Tests" begin
    include("atmotest.jl")
end

@testset "Rollout Tests" begin
    include("rollouttest.jl")
end
