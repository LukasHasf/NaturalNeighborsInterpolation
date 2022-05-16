include("../NaturalNeighborsInterpolation.jl")
using Plots

pointlist = randn(Float64, 2, 100)
pointlist[:, 1] = [2, -100]
pointlist[:, 2] = [1, -100]
scatter(pointlist)
savefig("pointlist.png")
display(pointlist)
hull = NaturalNeighborsInterpolation.ConvexHull2D(pointlist)