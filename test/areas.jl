using GeometricalPredicates
using Test
using Random
include("../NaturalNeighborsInterpolation.jl")

function dist(p, q)
    return sqrt((getx(p) - getx(q))^2 + (gety(p) - gety(q))^2)
end

function triArea(triangle)
    a,b,c = geta(triangle), getb(triangle), getc(triangle)
    ba = Point(getx(a) - getx(b), gety(a) - gety(b))
    bc = Point(getx(c) - getx(b), gety(c) - gety(b))
    area = 0.5 * (getx(ba) * gety(bc) - gety(ba) * getx(bc))
    return abs(area)
end

min_coord = 1.0 +eps(Float64)
max_coord = 2 - 2*eps(Float64)
width = max_coord - min_coord
n = 30
pointlist = min_coord .+ rand(Float64, 2, n) .* width
pointlist[:, 1] = [min_coord, min_coord]
pointlist[:, 2] = [min_coord, max_coord]
pointlist[:, 3] = [max_coord, max_coord]
points = Point2D[Point(pointlist[1,i], pointlist[2,i]) for i in 1:n]
triangles = [Primitive(points[i], points[i+1], points[i+2]) for i in 1:n-2]
for t in triangles
    correctArea = triArea(t)
    otherArea = NaturalNeighborsInterpolation.getArea([getpoint(t) for getpoint in [geta, getb, getc]])
    @test correctArea ≈ otherArea
end

square = [Point(min_coord, min_coord), Point(min_coord, max_coord),
            Point(max_coord, max_coord), Point(max_coord, min_coord)]
square = Random.shuffle(square)
correctArea = (max_coord - min_coord)^2
calculatedArea = NaturalNeighborsInterpolation.getArea(square)
@test correctArea ≈ calculatedArea
