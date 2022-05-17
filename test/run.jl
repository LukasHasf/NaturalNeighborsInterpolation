include("../NaturalNeighborsInterpolation.jl")
using VoronoiDelaunay
using Plots

min_coord = 1.0 +eps(Float64)
max_coord = 2 - 2*eps(Float64)
width = max_coord - min_coord

pointlist = min_coord .+ rand(Float64, 2, 64) .* width
pointlist[:, 1] .= [1.501, 1.5]
values = [pointlist[1,i] + pointlist[2,i] for i in 1:size(pointlist)[2]]
hull = NaturalNeighborsInterpolation.interpolate(pointlist, values, Point(1.5, 1.5))