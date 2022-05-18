include("../NaturalNeighborsInterpolation.jl")
using VoronoiDelaunay
using Plots
using ProgressBars

min_coord = 1.0 +eps(Float64)
max_coord = 2 - 2*eps(Float64)
width = max_coord - min_coord

#pointlist = min_coord .+ rand(Float64, 2, 128) .* width
pointlist = [1.45 1.6 1.4 1.6 1.4 1.6 1.55 1.7  1.4; 
             1.4  1.4 1.6 1.6 1.55 1.5 1.4  1.55 1.4] 
values = [pointlist[1,i] + pointlist[2,i] for i in 1:size(pointlist)[2]]
points = Point2D[Point(pointlist[1,i], pointlist[2,i]) for i in 1:size(pointlist)[2]]
#intp_value = NaturalNeighborsInterpolation.interpolate(pointlist, values, Point(1.5, 1.65))
# =3.2014440433213080
intp = NaturalNeighborsInterpolation.NaturalNeighborsInterpolator(points, values)
intp(Point(1.49, 1.56))

#= =#
xs = min_coord:0.01:max_coord
ys = min_coord:0.01:max_coord
#xs = 1.4:0.01:1.6
#ys = 1.4:0.01:1.6
interpolated_grid = Array{Float64, 2}(undef, length(ys), length(xs))

@simd for j in eachindex(xs)
@simd for i in eachindex(ys)
        print(xs[j], ", ", ys[i], "\r")
        interpolated_grid[i,j] = intp(Point(xs[j],ys[i]))
    end
end

heatmap(interpolated_grid)
savefig("interpolated_grid.png") # =#


