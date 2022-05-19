include("../src/NaturalNeighborsInterpolation.jl")
using VoronoiDelaunay
using Plots


#pointlist = min_coord .+ rand(Float64, 2, 128) .* width
pointlist = [0.45 0.6 0.4 0.6 0.4 0.6 0.55 0.7  0.4; 
             0.4  0.4 0.6 0.6 0.55 0.5 0.4  0.55 0.4] .+ 1
values = [pointlist[1,i] + pointlist[2,i] for i in 1:size(pointlist)[2]]
points = Point2D[Point(pointlist[1,i], pointlist[2,i]) for i in 1:size(pointlist)[2]]
intp = NaturalNeighborsInterpolation.NaturalNeighborsInterpolator(points, values)
xs = 0.3:0.01:0.7
ys = 0.3:0.01:0.7


#
xs = 1:0.01:2
ys = 1:0.01:2

interpolated_grid = Array{Float64, 2}(undef, length(ys), length(xs))

for j in eachindex(xs)
    for i in eachindex(ys)
        interpolated_grid[i,j] = intp(Point(xs[j],ys[i]))
    end
end

heatmap(interpolated_grid)
savefig("interpolated_grid.png") # =#


