using XEB
using Test

@testset "rewrite.jl" begin
    include("rewrite.jl")
end

@testset "Tensor networks" begin
    include("tensor_networks.jl")
end