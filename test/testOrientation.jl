using GeometricalPredicates
include("../NaturalNeighborsInterpolation.jl")


points = [Point(1.1, 1.1), Point(1.1, 1.3), Point(1.3, 1.3)][end:-1:1]
NaturalNeighborsInterpolation.polygonOrientation(poly)


poly = [Point(1.5, 1.75), Point(1.5, 1.25), Point(1.75, 1.5), Point(1.25, 1.5)]
poly_parts = NaturalNeighborsInterpolation.decomposeIntersection(poly)
for p in poly_parts
    println(NaturalNeighborsInterpolation.polygonOrientation(p))
end