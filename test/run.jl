include("../NaturalNeighborsInterpolation.jl")
using Plots

min_coord = 1.0 +eps(Float64)
max_coord = 2 - 2*eps(Float64)
width = max_coord - min_coord

pointlist = min_coord .+ rand(Float64, 2, 100) .* width
scatter(pointlist)
savefig("pointlist.png")
display(pointlist)
hull = NaturalNeighborsInterpolation.BowyerWatson(pointlist)