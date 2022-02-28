using Test
using CPEG
using BenchmarkTools
using StaticArrays
using LinearAlgebra
using SparseArrays

import CPEG

@testset "Atmosphere" begin
    include("atmotest.jl")
end

@testset "Gravity" begin
    include("gravtest.jl")
end

@testset "Rollout" begin
    include("rollouttest.jl")
end

@testset "Scaling" begin
    include("scalingtest.jl")
end

@testset "QP" begin
    include("qptest.jl")
end
