using NaturalNeighborsInterpolation
using VoronoiDelaunay
using Plots

#pointlist = min_coord .+ rand(Float64, 2, 128) .* width
pointlist = [0.45 0.6 0.4 0.6 0.4 0.6 0.55 0.7  0.4; 
             0.4  0.4 0.6 0.6 0.55 0.5 0.4  0.55 0.4] .+ 1
values = [pointlist[1,i] + pointlist[2,i] for i in 1:size(pointlist)[2]]
intp = NaturalNeighborsInterpolator(pointlist, values)

xs = 1.1:0.01:1.9
ys = 1.1:0.01:1.9

interpolated_grid = interpolateGrid(xs, ys, intp)
heatmap(interpolated_grid)
savefig("interpolated_grid.png")

intp_natural = NaturalNeighborsInterpolator(pointlist, values, fallback="natural")
interpolated_grid = interpolateGrid(xs, ys, intp_natural)
heatmap(interpolated_grid)
savefig("interpolated_grid_natural.png")

intp_natural = NaturalNeighborsInterpolator(pointlist, values, fallback="nan")
interpolated_grid = interpolateGrid(xs, ys, intp_natural)
heatmap(interpolated_grid)
savefig("interpolated_grid_nan.png")
