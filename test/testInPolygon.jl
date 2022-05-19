using GeometricalPredicates
using Plots
include("../src/NaturalNeighborsInterpolation.jl")

xs = 1.1:0.01:1.9
ys = 1.1:0.01:1.9

cell = [  Point2D(1.5, 1.5750000000000004)
Point2D(1.5, 1.5749999999999997)
Point2D(1.4550000000000003, 1.4849999999999997)
Point2D(1.5, 1.4343750000000004)
Point2D(1.5749999999999997, 1.4437499999999999)
Point2D(1.5833333333333337, 1.4500000000000002)
Point2D(1.5500000000000005, 1.5499999999999998)]
cell = unique(cell)


function reduce(cell)
reduced_cell = []
checked_points = []
for p1 in cell
    isOkay = true
    for p2 in cell
        if p1==p2
            continue
        end
        if NaturalNeighborsInterpolation.distSq(p1, p2) < 1e-30 && !(p2 in checked_points)
            println("Wrong")
            isOkay = false
        end
    end
    if isOkay
        push!(reduced_cell, p1)
    end
    push!(checked_points, p1)
end
return reduced_cell
end



xq = [getx(p) for p in cell]
yq = [gety(p) for p in cell]
points = [[x,y] for (x,y) in zip(xq,yq)]
display(unique(points))


polygon = Polygon(cell...)

mask1 = zeros(Int, length(ys), length(xs))
mask2 = zeros(Int, length(ys), length(xs))
for (j,x) in enumerate(xs)
    for (i,y) in enumerate(ys)
        mask1[i,j] = inpolygon(polygon, Point(x,y)) ? 1 : 0
        mask2[i,j] = NaturalNeighborsInterpolation.inPolygon(cell, Point(x,y)) ? 1 : 0
    end
end

heatmap(mask1)
savefig("inpolygon.svg")

heatmap(mask2)
savefig("inPolygon.svg")