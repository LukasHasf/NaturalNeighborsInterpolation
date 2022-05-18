include("../src/NaturalNeighborsInterpolation.jl")
using Test

@testset "Check area calculation" begin
    include("areas.jl")
end

@testset "Check self intersection detection" begin
    include("lines.jl")
end