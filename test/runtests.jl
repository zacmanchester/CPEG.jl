using Test
using CPEG
using BenchmarkTools
using StaticArrays
using LinearAlgebra
using SparseArrays

import CPEG

@testset "Atmosphere Tests" begin
    include("atmotest.jl")
end

@testset "Rollout Tests" begin
    include("rollouttest.jl")
end

@testset "Scaling Tests" begin
    include("scalingtest.jl")
end

@testset "QP Tests" begin
    include("qptest.jl")
end
