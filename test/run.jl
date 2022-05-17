include("../NaturalNeighborsInterpolation.jl")
using VoronoiDelaunay
using Plots

min_coord = 1.0 +eps(Float64)
max_coord = 2 - 2*eps(Float64)
width = max_coord - min_coord

#pointlist = min_coord .+ rand(Float64, 2, 64) .* width
pointlist = [1.45 1.6 1.4 1.6 1.4 1.6 1.55 1.7  1.4; 
             1.4  1.4 1.6 1.6 1.55 1.5 1.4  1.55 1.4] 
#pointlist .+= rand(Float64, size(pointlist)) * 0.001
values = [pointlist[1,i] + pointlist[2,i] for i in 1:size(pointlist)[2]]
hull = NaturalNeighborsInterpolation.interpolate(pointlist, values, Point(1.5, 1.5))